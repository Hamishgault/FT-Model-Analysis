param(
    [string]$Python = "python"
)

Write-Host "Creating virtual environment at .venv..."
& $Python -m venv .venv

Write-Host "Upgrading pip and installing requirements..."
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\pip.exe install -r requirements.txt

Write-Host "Done. To activate the environment in PowerShell run: . .venv\\Scripts\\Activate.ps1"
