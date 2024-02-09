sphinx-apidoc -f -o docs tivatci
call %~dp0\docs\make.bat html
start chrome "%~dp0/docs/index.html"
PAUSE
