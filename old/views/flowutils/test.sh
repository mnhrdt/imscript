./viewflow -1 hom.flo hom_colors.png
./flowarrows 0.1 19 hom.flo hom_arrows.png
./backflow hom.flo lena.png hom_lena.png
./flowinv 10 0.01 hom.flo ihom.flo
./viewflow -1 ihom.flo hom_colors.png
./flowarrows 0.1 19 ihom.flo hom_arrows.png
./backflow hom.flo lena.png ihom_lena.png
