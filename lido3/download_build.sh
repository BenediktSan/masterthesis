#!/bin/bash -l

scp -i  ~/.ssh/id_rsa -r smbtsand@gw01.lido.tu-dortmund.de:/work/smbtsand/Masterarbeit/main/build ~/Data/Masterarbeit/main
scp -i  ~/.ssh/id_rsa -r smbtsand@gw01.lido.tu-dortmund.de:/work/smbtsand/Masterarbeit/lido3 ~/Data/Masterarbeit
