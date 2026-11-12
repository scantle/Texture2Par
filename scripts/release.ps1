[CmdletBinding()]
param(
    [Parameter(Mandatory = $true, Position = 0)]
    [ValidatePattern('^v?\d+\.\d+\.\d+(?:[-+][0-9A-Za-z.-]+)?$')]
    [string]$Version,

    [switch]$SkipTests,
    [switch]$Publish,
    [switch]$Draft,
    [switch]$Prerelease
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Program,

        [Parameter(ValueFromRemainingArguments = $true)]
        [string[]]$Arguments
    )

    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Program $($Arguments -join ' ')"
    }
}

function Get-DevenvPath {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (Test-Path -LiteralPath $vswhere) {
        $installation = & $vswhere -latest -products '*' -property installationPath |
            Select-Object -First 1
        if ($installation) {
            $candidate = Join-Path $installation 'Common7\IDE\devenv.com'
            if (Test-Path -LiteralPath $candidate) {
                return $candidate
            }
        }
    }

    $knownEditions = @('Community', 'Professional', 'Enterprise')
    foreach ($edition in $knownEditions) {
        $candidate = Join-Path ${env:ProgramFiles} "Microsoft Visual Studio\2022\$edition\Common7\IDE\devenv.com"
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
    }

    throw 'Visual Studio 2022 devenv.com was not found.'
}

function Get-OpenMPRuntimePath {
    if ($env:CMPLR_ROOT) {
        $candidate = Join-Path $env:CMPLR_ROOT 'bin\libiomp5md.dll'
        if (Test-Path -LiteralPath $candidate) {
            return $candidate
        }
    }

    $roots = @(
        (Join-Path ${env:ProgramFiles(x86)} 'Intel\oneAPI\compiler'),
        (Join-Path ${env:ProgramFiles} 'Intel\oneAPI\compiler')
    )
    $candidate = Get-ChildItem -LiteralPath $roots -Filter 'libiomp5md.dll' -Recurse -File -ErrorAction SilentlyContinue |
        Where-Object { $_.DirectoryName -match '[\\/]bin$' } |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if ($candidate) {
        return $candidate.FullName
    }

    throw 'Intel OpenMP runtime libiomp5md.dll was not found.'
}

function Assert-PathWithin {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path,

        [Parameter(Mandatory = $true)]
        [string]$Root
    )

    $fullPath = [IO.Path]::GetFullPath($Path)
    $fullRoot = [IO.Path]::GetFullPath($Root).TrimEnd([IO.Path]::DirectorySeparatorChar)
    if (-not $fullPath.StartsWith($fullRoot + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to modify path outside release directory: $fullPath"
    }
}

$normalizedVersion = $Version.TrimStart('v')
$tag = "v$normalizedVersion"
$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$solution = Join-Path $repoRoot 'Texture2Par.sln'
$versionFile = Join-Path $repoRoot 'version.f90'
$binary = Join-Path $repoRoot 'bin\Texture2Par.exe'
$binaryRuntime = Join-Path $repoRoot 'bin\libiomp5md.dll'
$releaseRoot = Join-Path $repoRoot 'build\release'
$packageName = "Texture2Par-v$normalizedVersion-windows-x64"
$packageDirectory = Join-Path $releaseRoot $packageName
$archive = Join-Path $releaseRoot "$packageName.zip"
$checksum = "$archive.sha256"
$pytestTemp = Join-Path $releaseRoot 'pytest-tmp'

Set-Location -LiteralPath $repoRoot

if (($Draft -or $Prerelease) -and -not $Publish) {
    throw '-Draft and -Prerelease require -Publish.'
}
if ($Publish -and $SkipTests) {
    throw 'Published releases must run the test suite; remove -SkipTests.'
}

$versionText = Get-Content -LiteralPath $versionFile -Raw
$versionMatch = [regex]::Match($versionText, 'VERSION\s*=\s*"([^"]+)"')
if (-not $versionMatch.Success) {
    throw 'Could not read VERSION from version.f90.'
}
if ($versionMatch.Groups[1].Value -ne $normalizedVersion) {
    throw "Requested version $normalizedVersion does not match version.f90 ($($versionMatch.Groups[1].Value))."
}

$ghCommand = $null
if ($Publish) {
    $ghCommand = Get-Command gh -ErrorAction SilentlyContinue
    if (-not $ghCommand) {
        throw 'GitHub CLI (gh) is required for -Publish. Install it, authenticate with gh auth login, and rerun.'
    }

    $changes = @(& git status --porcelain)
    if ($LASTEXITCODE -ne 0) {
        throw 'Unable to inspect Git status.'
    }
    if ($changes.Count -ne 0) {
        throw 'Publishing requires a clean worktree. Commit or stash all changes first.'
    }
}

$devenv = Get-DevenvPath
Write-Host "Building $packageName..."
Invoke-Checked $devenv $solution '/Rebuild' 'Release_OpenMP|x64'

if (-not (Test-Path -LiteralPath $binary)) {
    throw "Expected build output was not created: $binary"
}

$runtime = Get-OpenMPRuntimePath
Copy-Item -LiteralPath $runtime -Destination $binaryRuntime -Force
New-Item -ItemType Directory -Path $releaseRoot -Force | Out-Null

if (-not $SkipTests) {
    Assert-PathWithin -Path $pytestTemp -Root $releaseRoot
    if (Test-Path -LiteralPath $pytestTemp) {
        Remove-Item -LiteralPath $pytestTemp -Recurse -Force
    }

    Write-Host "Testing $binary..."
    $previousT2PExe = $env:T2P_EXE
    try {
        $env:T2P_EXE = $binary
        $venvPython = Join-Path $repoRoot '.venv-tests\Scripts\python.exe'
        $conda = Get-Command conda -ErrorAction SilentlyContinue
        if (Test-Path -LiteralPath $venvPython) {
            Invoke-Checked $venvPython '-m' 'pytest' '-q' '-p' 'no:cacheprovider' "--basetemp=$pytestTemp" '-x' '--tb=long'
        }
        elseif ($conda) {
            Invoke-Checked $conda.Source 'run' '-n' 't2p_tests' 'python' '-m' 'pytest' '-q' '-p' 'no:cacheprovider' "--basetemp=$pytestTemp" '-x' '--tb=long'
        }
        else {
            throw 'Neither the t2p_tests Conda environment nor .venv-tests was found.'
        }
    }
    finally {
        $env:T2P_EXE = $previousT2PExe
    }
}
else {
    Write-Warning 'Tests were skipped.'
}

Write-Host "Packaging $archive..."
foreach ($path in @($packageDirectory, $archive, $checksum)) {
    Assert-PathWithin -Path $path -Root $releaseRoot
    if (Test-Path -LiteralPath $path) {
        Remove-Item -LiteralPath $path -Recurse -Force
    }
}

New-Item -ItemType Directory -Path $packageDirectory -Force | Out-Null
Copy-Item -LiteralPath $binary -Destination (Join-Path $packageDirectory 'Texture2Par.exe')
Copy-Item -LiteralPath $binaryRuntime -Destination (Join-Path $packageDirectory 'libiomp5md.dll')
Copy-Item -LiteralPath (Join-Path $repoRoot 'README.md') -Destination (Join-Path $packageDirectory 'README.md')

$packageNotes = @"
Texture2Par v$normalizedVersion for 64-bit Windows

ABOUT libiomp5md.dll

libiomp5md.dll is Intel's OpenMP runtime library. Texture2Par uses it to create,
schedule, and synchronize worker threads during parallel deterministic
libkriging solves. It is a required program dependency, not an input file or
optional plugin.

Keep Texture2Par.exe and libiomp5md.dll in the same directory. Windows searches
the executable's directory for dependent DLLs, including when Texture2Par is
launched from a different working directory. If you copy Texture2Par.exe into a
model or working directory, copy libiomp5md.dll with it.

The DLL is required even when NTHREADS is 1 because Windows loads it when the
program starts. If Windows reports that libiomp5md.dll is missing, re-extract
this complete ZIP or copy this bundled DLL beside Texture2Par.exe. Do not
download a replacement from a third-party DLL site.

NTHREADS defaults to 1. Set NTHREADS greater than 1 in the input file to enable
parallel deterministic libkriging solves.
"@
Set-Content -LiteralPath (Join-Path $packageDirectory 'DISTRIBUTION.txt') -Value $packageNotes -Encoding ASCII

Compress-Archive -LiteralPath $packageDirectory -DestinationPath $archive -CompressionLevel Optimal
$sha256 = [Security.Cryptography.SHA256]::Create()
$archiveStream = [IO.File]::OpenRead($archive)
try {
    $archiveHash = ([BitConverter]::ToString($sha256.ComputeHash($archiveStream))).Replace('-', '').ToLowerInvariant()
}
finally {
    $archiveStream.Dispose()
    $sha256.Dispose()
}
"$archiveHash  $([IO.Path]::GetFileName($archive))" |
    Set-Content -LiteralPath $checksum -Encoding ASCII

Write-Host "Created $archive"
Write-Host "Created $checksum"

if ($Publish) {
    $head = (& git rev-parse HEAD).Trim()
    if ($LASTEXITCODE -ne 0) {
        throw 'Unable to resolve Git HEAD.'
    }
    $existingTag = @(& git tag --list $tag)
    if ($LASTEXITCODE -ne 0) {
        throw "Unable to inspect Git tag $tag."
    }
    if ($existingTag.Count -eq 0) {
        Invoke-Checked 'git' 'tag' '-a' $tag '-m' "Texture2Par $normalizedVersion"
    }
    else {
        $tagCommit = (& git rev-list '-n' '1' $tag).Trim()
        if ($LASTEXITCODE -ne 0 -or $tagCommit -ne $head) {
            throw "Tag $tag already exists and does not point to HEAD."
        }
    }

    Invoke-Checked 'git' 'push' 'origin' 'HEAD'
    Invoke-Checked 'git' 'push' 'origin' $tag

    $releaseArguments = @(
        'release', 'create', $tag,
        $archive, $checksum,
        '--verify-tag',
        '--title', "Texture2Par $normalizedVersion",
        '--generate-notes'
    )
    if ($Draft) {
        $releaseArguments += '--draft'
    }
    if ($Prerelease) {
        $releaseArguments += '--prerelease'
    }
    Invoke-Checked $ghCommand.Source @releaseArguments
}

Write-Host 'Release workflow completed successfully.'
