rm -rf picture
mkdir picture

for i in 1 2 3 4 5 6 7 8 9
do
  scp yamad@venus:/home/yamad/HaploScore/figure_cut$i.png ./picture
  scp yamad@venus:/home/yamad/HaploScore/figure_shf$i.png ./picture
done
