@echo off
if not exist results mkdir results

bin\Density1DRods_VS.exe 2 40 5000 

echo Simulations started.
pause