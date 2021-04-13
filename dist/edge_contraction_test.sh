cd ..
make
cd ./dist
echo $1+"_"+$3 > ../result/edge_IA.txt
./MM_gradient -i /mnt/c/things/projects/Terrain_Trees/data/$1 -c $2 $3 >>  ../result/edge_IA.txt 2>&1
