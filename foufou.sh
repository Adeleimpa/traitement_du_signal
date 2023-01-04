echo "Script d'execution des divers cas Fourier"
# All 5 sec 44100Hz
# Gamme chroma en La sans filtre
echo Gamme Chroma No filtre
./fourier C 5 44100 0 0 0 > Chroma_5SEC_44100Hz_Filter_none.log
# Gamme chroma filtre bas, haut , 600Hz butter 450Hz
echo Gamme Chroma Filtre passe-bas 600Hz
./fourier C 5 44100 1 600 0 > Chroma_5SEC_44100Hz_Filter_Low_600.log
echo Gamme Chroma Filtre passe-haut 600Hz
./fourier C 5 44100 2 600 0 > Chroma_5SEC_44100Hz_Filter_High_600.log
echo Gamme Chroma Filtre Butter 450Hz
./fourier C 5 44100 3 450 0 > Chroma_5SEC_44100Hz_Filter_Butter_450.log
echo Gamme Chroma Filtre Coupe Bande 500Hz - 550Hz
./fourier C 5 44100 5 500 550 > Chroma_5SEC_44100Hz_Filter_CutBand_500_550.log
# Accord domisol filtre bas, haut 350Hz, butter 300Hz
echo Accord Filtre passe-bas 350Hz
./fourier A 5 44100 1 350 0 > Accord_5SEC_44100Hz_Filter_Low_350.log
echo Accord Filtre passe-haut 350Hz
./fourier A 5 44100 2 350 0 > Accord_5SEC_44100Hz_Filter_High_350.log
echo Accord Filtre Butter 300Hz
./fourier A 5 44100 3 300 0 > Accord_5SEC_44100Hz_Filter_Butter_300.log
# Accord C2 generated
echo Accord generated
./fourier C2 5 44100 0 0 0 > Accord_GEN_5SEC_44100Hz.log
# Gamme Piano Bas, Haut 350Hz, Butter 350Hz, Pass-Band 300-400
echo Piano Filtre passe-bas 350Hz
./fourier P 5 44100 1 350 0 > Piano_Filter_Low_350.log
echo Piano Filtre passe-haut 350Hz
./fourier P 5 44100 2 350 0 > Piano_Filter_High_350.log
echo Piano Filtre Butter 350Hz
./fourier P 5 44100 3 350 0 > Piano_Filter_Butter_350.log
echo Piano Filter Pass-Band 300-400
./fourier P 5 44100 4 300 400 > Piano_Filter_Pass_Band_300_400.log
echo Piano Filter Cut-Band 300-400
./fourier P 5 44100 5 300 400 > Piano_Filter_Cut_Band_300_400.log
