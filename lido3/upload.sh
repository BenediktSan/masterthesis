#!/bin/bash -l

scp -i  ~/.ssh/id_rsa -r ~/Data/Masterarbeit/lido3 smbtsand@gw01.lido.tu-dortmund.de:/work/smbtsand/Masterarbeit
scp -i  ~/.ssh/id_rsa -r ~/Data/Masterarbeit/header smbtsand@gw01.lido.tu-dortmund.de:/work/smbtsand/Masterarbeit
scp -i  ~/.ssh/id_rsa -r ~/Data/Masterarbeit/main smbtsand@gw01.lido.tu-dortmund.de:/work/smbtsand/Masterarbeit
