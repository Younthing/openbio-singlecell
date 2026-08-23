[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python,
    [string]$ComfyRoot,
    [string]$FrontendRoot,
    [Alias("h")][switch]$Help,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$ComfyArgs
)

$ErrorActionPreference = "Stop"

if ($Help) {
    Write-Output @"
Usage: powershell -File scripts/start.ps1 [-Python <python.exe>] [-ComfyRoot <path>] [-FrontendRoot <path>] [ComfyUI arguments...]

Starts ComfyUI with the sibling OpenBio frontend dist, local assets,
and API nodes disabled. Additional arguments are forwarded to ComfyUI.

ComfyUI can also be selected with OPENBIO_COMFYUI_ROOT. The frontend repository
can also be selected with OPENBIO_FRONTEND_ROOT.
"@
    exit 0
}

$PluginRoot = Split-Path -Parent $PSScriptRoot

function Resolve-OpenBioComfyRoot([string]$Requested) {
    if ($Requested) {
        $Candidates = @($Requested)
        $Source = "-ComfyRoot"
    } elseif ($env:OPENBIO_COMFYUI_ROOT) {
        $Candidates = @($env:OPENBIO_COMFYUI_ROOT)
        $Source = "OPENBIO_COMFYUI_ROOT"
    } else {
        $Candidates = @()
        $PluginParent = Split-Path -Parent $PluginRoot
        if ((Split-Path -Leaf $PluginParent) -eq "custom_nodes") {
            $Candidates += Split-Path -Parent $PluginParent
        }
        $Candidates += Join-Path $PluginParent "ComfyUI"
        $Source = $null
    }

    foreach ($Candidate in $Candidates) {
        if (Test-Path -LiteralPath $Candidate -PathType Container) {
            $Resolved = (Resolve-Path -LiteralPath $Candidate).Path
            if (Test-Path -LiteralPath (Join-Path $Resolved "main.py") -PathType Leaf) {
                return $Resolved
            }
        }
    }

    if ($Source) {
        throw "$Source does not point to a ComfyUI root containing main.py: $($Candidates[0])"
    }
    throw "ComfyUI root was not found. Pass -ComfyRoot or set OPENBIO_COMFYUI_ROOT."
}

function Resolve-OpenBioFrontendRoot([string]$Requested) {
    if ($Requested) {
        $Candidates = @($Requested)
        $Source = "-FrontendRoot"
    } elseif ($env:OPENBIO_FRONTEND_ROOT) {
        $Candidates = @($env:OPENBIO_FRONTEND_ROOT)
        $Source = "OPENBIO_FRONTEND_ROOT"
    } else {
        $Candidates = @(
            (Join-Path (Split-Path -Parent $PluginRoot) "ComfyUI_frontend"),
            (Join-Path (Split-Path -Parent $ComfyRoot) "ComfyUI_frontend")
        )
        $Source = $null
    }

    foreach ($Candidate in $Candidates) {
        if (Test-Path -LiteralPath $Candidate -PathType Container) {
            return (Resolve-Path -LiteralPath $Candidate).Path
        }
    }

    if ($Source) {
        throw "$Source does not point to a frontend repository: $($Candidates[0])"
    }
    throw "ComfyUI frontend repository was not found. Pass -FrontendRoot or set OPENBIO_FRONTEND_ROOT."
}

$ComfyRoot = Resolve-OpenBioComfyRoot $ComfyRoot
$FrontendRoot = Resolve-OpenBioFrontendRoot $FrontendRoot
$FrontendDist = Join-Path $FrontendRoot "dist"

$RequiredFrontendFiles = @("index.html", "LICENSE", "THIRD_PARTY_NOTICES.md")
$MissingFrontendFiles = @($RequiredFrontendFiles | Where-Object {
    $Candidate = Join-Path $FrontendDist $_
    (-not (Test-Path -LiteralPath $Candidate -PathType Leaf)) -or ((Get-Item -LiteralPath $Candidate).Length -eq 0)
})
if ($MissingFrontendFiles.Count -gt 0) {
    throw @"
OpenBio frontend dist is incomplete: $FrontendDist
Missing: $($MissingFrontendFiles -join ", ")
Build it first:
  Set-Location "$FrontendRoot"
  corepack pnpm install --frozen-lockfile
  corepack pnpm build:openbio
  corepack pnpm exec node "$PluginRoot/scripts/generate_frontend_notices.mjs"
"@
}

function Resolve-OpenBioPython([string]$Requested) {
    $Candidates = if ($Requested) {
        @($Requested)
    } else {
        @(
            (Join-Path $ComfyRoot ".venv\Scripts\python.exe"),
            (Join-Path $ComfyRoot "venv\Scripts\python.exe"),
            "python"
        )
    }

    foreach ($Candidate in $Candidates) {
        if (Test-Path -LiteralPath $Candidate) {
            return (Resolve-Path -LiteralPath $Candidate).Path
        }
        $Command = Get-Command $Candidate -ErrorAction SilentlyContinue
        if ($Command) {
            return $Command.Source
        }
    }
    throw "Python was not found. Install Python 3.12+ or pass -Python with the ComfyUI interpreter."
}

$PythonExe = Resolve-OpenBioPython $Python
& $PythonExe -c "import sys; raise SystemExit(sys.version_info < (3, 12))"
if ($LASTEXITCODE -ne 0) {
    throw "openbio-singlecell requires Python 3.12 or newer."
}

$Arguments = @(
    (Join-Path $ComfyRoot "main.py"),
    "--disable-api-nodes",
    "--enable-assets",
    "--front-end-root",
    $FrontendDist
)
$Arguments += $ComfyArgs

Write-Host "Using Python: $PythonExe"
Push-Location $ComfyRoot
try {
    & $PythonExe @Arguments
    $ExitCode = $LASTEXITCODE
} finally {
    Pop-Location
}
exit $ExitCode
