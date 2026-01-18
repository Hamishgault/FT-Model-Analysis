param(
    [string]$Python = ".venv\Scripts\python.exe"
)

Write-Host "Activating venv and running pytest..."
. .venv\Scripts\Activate.ps1
pytest --cov=./ -q
