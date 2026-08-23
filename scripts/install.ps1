param(
    [string]$Python,
    [string]$ComfyRoot,
    [switch]$ForceDemo,
    [Alias("h")][switch]$Help
)

$ErrorActionPreference = "Stop"

if ($Help) {
    Write-Output @"
Usage: powershell -File scripts/install.ps1 [-Python <python.exe>] [-ComfyRoot <path>] [-ForceDemo]

Installs only openbio-singlecell Python requirements, checks imports, and
generates the local demonstration AnnData H5AD file.

ComfyUI is resolved from -ComfyRoot, OPENBIO_COMFYUI_ROOT, the Manager install
location, or a sibling ComfyUI repository (in that order).
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

$ComfyRoot = Resolve-OpenBioComfyRoot $ComfyRoot

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

function Invoke-CheckedPython([string[]]$Arguments) {
    & $PythonExe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Python command failed with exit code $LASTEXITCODE."
    }
}

$PythonExe = Resolve-OpenBioPython $Python
Write-Host "Using Python: $PythonExe"
& $PythonExe -c "import sys; raise SystemExit(sys.version_info < (3, 12))"
if ($LASTEXITCODE -ne 0) {
    throw "openbio-singlecell requires Python 3.12 or newer."
}
$UvCommand = Get-Command "uv" -ErrorAction SilentlyContinue
if ($UvCommand) {
    & $UvCommand.Source pip install --python $PythonExe -r (Join-Path $PluginRoot "requirements.txt")
    if ($LASTEXITCODE -ne 0) {
        throw "uv dependency installation failed with exit code $LASTEXITCODE."
    }
} else {
    & $PythonExe -m pip --version *> $null
    if ($LASTEXITCODE -ne 0) {
        throw "Neither uv nor pip is available. Install uv or add pip to the ComfyUI Python environment."
    }
    Invoke-CheckedPython @("-m", "pip", "install", "--disable-pip-version-check", "-r", (Join-Path $PluginRoot "requirements.txt"))
}

Invoke-CheckedPython @((Join-Path $PSScriptRoot "check_imports.py"), $PluginRoot, $ComfyRoot)

$DemoArguments = @((Join-Path $PSScriptRoot "generate_demo.py"), "--comfy-root", $ComfyRoot)
if ($ForceDemo) {
    $DemoArguments += "--force"
}
Invoke-CheckedPython $DemoArguments
Write-Host "openbio-singlecell installation complete."
