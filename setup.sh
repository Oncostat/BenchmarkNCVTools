
mkdir databases results scores



g++ -std=c++11 src/coordfusion.cpp -o src/coordfusion
g++ -std=c++11 src/snpsomDM.cpp -o src/snpsomDM
g++ -std=c++11 src/cosmicDM.cpp -o src/cosmicDM
g++ -std=c++11 src/rmUseless.cpp -o src/rmUseless
g++ -std=c++11 src/rmMultiDupl.cpp -o src/rmMultiDupl
g++ -std=c++11 src/clinvarDM.cpp -o src/clinvarDM
g++ -std=c++11 src/exportGWAVA.cpp -o src/exportGWAVA
g++ -std=c++11 src/getGWAVAsc.cpp -o src/getGWAVAsc
g++ -std=c++11 src/makeIndex.cpp -o src/makeIndex
g++ -std=c++11 src/variantScoring.cpp -o src/variantScoring
g++ -std=c++11 src/mergeSnpsom.cpp -o src/mergeSnpsom
g++ -std=c++11 src/filterRec.cpp -o src/filterRec
g++ -std=c++11 src/filterAF.cpp -o src/filterAF
g++ -std=c++11 src/sampling.cpp -o src/sampling
g++ -std=c++11 src/region.cpp -o src/region
g++ -std=c++11 src/uniquePos.cpp -o src/uniquePos
