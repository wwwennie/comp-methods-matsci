import numpy as np

def fccmke(nc):
   """ Generate fractional coordinates of atoms in FCC supercell
       Input:
           nc (integer): supercell of nc x nc x nc dimensions
       Output:
           s (array): fractional atomic coordinates of supercell
   """ 

   natoms = 4
   r = np.array([[0, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
   i1 = 0
   s = np.zeros((natoms * nc**3, 3))

   # in fractional coordinates
   for k in range(1, nc + 1):
      for l in range(1, nc + 1):
         for m in range(1, nc + 1):
            for i in range(natoms):
               s[i1, 0] = (r[i, 0] + k - 1) / nc
               s[i1, 1] = (r[i, 1] + l - 1) / nc
               s[i1, 2] = (r[i, 2] + m - 1) / nc
               i1 += 1
   return s

