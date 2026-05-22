# A simple batch file to call gforran and compile L_RANGE

del *.o

"C:\gcc\bin\gfortran" -c  Parm_Vars_and_Vals.f95
"C:\gcc\bin\gfortran" -c  Structures.f95
"C:\gcc\bin\gfortran" -c  State_Variables.f95
"C:\gcc\bin\gfortran" -c  Decomposition.f95
"C:\gcc\bin\gfortran" -c  Initialize_Model.f95
"C:\gcc\bin\gfortran" -c  Misc_Material.f95
"C:\gcc\bin\gfortran" -c  Outputs.f95
"C:\gcc\bin\gfortran" -c  Plant_Death.f95
"C:\gcc\bin\gfortran" -c  Plant_Populations.f95
"C:\gcc\bin\gfortran" -c  Productivity.f95
"C:\gcc\bin\gfortran" -c  Soil_and_Water.f95
"C:\gcc\bin\gfortran" -c  Weather.f95
"C:\gcc\bin\gfortran" -c  Herbivore_Populations.f95

"C:\gcc\bin\gfortran" -static-libgfortran -static-libgcc Af_Range.f95 Decomposition.o Initialize_Model.o Misc_Material.o Outputs.o Parm_Vars_and_Vals.o Plant_Death.o Plant_Populations.o Productivity.o Soil_and_Water.o State_Variables.o Structures.o Weather.o Herbivore_Populations.o -o F:\AF_Range\AF_Range_Bin\Af_Range.exe

