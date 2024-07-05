library(bio3d)
dcd <- read.dcd("D:/Abel/Nitte_colab/PCA/PCA_DCCM/MD_center_neratinib_1000.dcd")
pdb <- read.pdb("D:/Abel/Nitte_colab/PCA/MD_neratinib.pdb")
ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)
dim(xyz) == dim(dcd)
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )
cij<-dccm(xyz[,ca.inds$xyz])
plot(cij)
pymol(cij, pdb, type="launch",exefile="C:/ProgramData/pymol/PyMOLWin.exe")
