cd ..
make
cd ./dist
echo $1+"_"+$2 > ../result/experiments_IA.txt
./MM_gradient -i /mnt/d/things/projects/Terrain_Trees/data/$1 -s $2 >>  ../result/experiments_IA.txt 2>&1
