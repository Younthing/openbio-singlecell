function Resolve-OpenBioComfyRoot(
    [string]$Requested,
    [string]$PluginRoot
) {
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

function Resolve-OpenBioPython(
    [string]$Requested,
    [string]$ComfyRoot
) {
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
