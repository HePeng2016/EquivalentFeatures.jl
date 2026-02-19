# Introduction
To do 

     using  EquivalentFeatures 
     
# Initial
     EquivalentFeatures.setN(N_new);
*N_new* indicates the maximum orbital quantum number *l* of the input spherical harmonic representation *Y(lm)*.

Because in julia the index begins from 1, so we define *N_new* as *l+1*. In the following, the definition of degree n is the quantum number *l* plus one. If this command is missing, then the default maximum orbital quantum number *l* is 2, so *N_new* is 3. 
     
     EquivalentFeatures.Initial();
   Once this function had been carried out, the tables for the Clebsch Gordan coefficients and the Wigner 3j coefficients were generated. And the data structures for storing the relations between coefficients and variables were generated.
# Usage
    EquivalentFeatures.CtoS_Encode(V1);

  This function changes the cartesian coordination terms to the spherical harmonic tensor terms. 
  
  *V1* is the vector that records the multipole moment terms of the multipole expansion in the Cartesian system of coordinates.
  The order of the terms in this vector is determined by the length of the term, if the terms are equally long, they are ordered according to the lexicographical order. 
  
  e.g.
  
   *x<y<z*,
   *xx<x* ,
   *xx<xy<xz<yy<yz<zz*
   
*V1* is formatted as: 

N=1: [monopole]

*[1]*
    
N=2: [monopole,dipole]

*[1,x,y,z]*
    
N=3: [monopole,dipole,quadruple]

*[1,x,y,z,xx,xy,xz,yy,yz,zz]*

The output is a spherical harmonic tensor stored as a complex vector.

e.g.

N=3:

*[Y(l=0,m=0),Y(l=1,m=-1),Y(l=1,m=0),Y(l=1,m=1),Y(l=2,m=-2),Y(l=2,m=-1),Y(l=2,m=0),Y(l=2,m=1),Y(l=2,m=2)]*

*Y(l,m)* is the term for spherical multipole moment in the spherical coordinate system. *l* is the orbital quantum number and *m* is the azimuthal quantum number.

     EquivalentFeatures.IncreaseDegree(S2,n);
     
This function increases the degree of the spherical harmonic tensor by n. S2 should be a spherical harmonic tensor [monopole,dipole] derived from a normalized Cartesian coordinate vector.

     EquivalentFeatures.SelfProduct(V1)

This function will change a spherical harmonic tensor to a rotation invariance embedding vector.

     EquivalentFeatures.SelfProduct(V1,n,n2)

For the Clebsch–Gordan (CG) product of any two shells in the V1 spherical harmonic tensor, only top n2  with smaller degree are selected and the maximum degree is n.  
     
     EquivalentFeatures.SelfProductPairwise(V1,n,n2)
     
The definition of n,n2 is the same as SelfProduct function. For the rotation-invariance embedding, both norm and pairwise product are considered. 

Encoding one object with the same center:

e.g. 

      using  EquivalentFeatures
      using LinearAlgebra
	  
      EquivalentFeatures.Initial();  
      C1 =  [0.0043477849927746155,0.0,0.9999905483381614];
      R1 =  norm(C1,2);
      C1 =  C1/R1; # Cartesian coordinate normalization.
      V1 =  vcat(1, C1);
      
      C2 =  [0.772027518982468,0.33454525822573616,0.5404192632877276];
      R2 =  norm(C2,2);
      C2 =  C2/R2;# Cartesian coordinate normalization.
      V2 =  vcat(1, C2);
      # c2 is the vector obtained by rotating c1.
	  # V1 and V2 are two identical spherical harmonic tensors with different rotations. 
      
      S1 = EquivalentFeatures.CtoS_Encode(V1,2);
      S2 = EquivalentFeatures.CtoS_Encode(V2,2);    # Convert the vector to a dipole spherical harmonic tensor(*l=1*), because n indicates *l+1*, the input *n = 2*.  
      S1 = EquivalentFeatures.IncreaseDegree(S1,1);  
      S2 = EquivalentFeatures.IncreaseDegree(S2,1); # Increase the dipole spherical harmonic tensor into quadrupole spherical harmonic tensor. 
	    S1[2:4] = S1[2:4]*R1; 
	    S1[5:9] = S1[5:9]*R1*R1; 
	    S2[2:4] = S2[2:4]*R2; 
	    S2[5:9] = S2[5:9]*R2*R2; 
      E1 = EquivalentFeatures.SelfProduct(S1); 
      E2 = EquivalentFeatures.SelfProduct(S2);
      Loss = sum(abs.(E1 - E2))
      E1 = EquivalentFeatures.SelfProduct(S1,2,2); 
      E2 = EquivalentFeatures.SelfProduct(S2,2,2);
      Loss = sum(abs.(E1 - E2))
      E1 = EquivalentFeatures.SelfProductPairwise(S1,2,2);
      E2 = EquivalentFeatures.SelfProductPairwise(S2,2,2);
      Loss = sum(abs.(E1 - E2)) 



Encoding two objects with the same center: 


      EquivalentFeatures.ProductEncode(V1,V2)


This function will return the invariant coding of the tensor product of two spherical harmonic tensors (V1, V2).

e.g.

      V1 =[1.0,0.0043477849927746155,0.0,0.9999905483381614,0.11]; 
      V2 =[1.0,0.772027518982468,0.33454525822573616,0.5404192632877276];

      using  EquivalentFeatures
      EquivalentFeatures.setN(2); 
      EquivalentFeatures.Initial();
      S1 = EquivalentFeatures.CtoS_Encode(V1);
      S2 = EquivalentFeatures.CtoS_Encode(V2);
      EquivalentFeatures.ProductEncode(S1,S2);
S1 and S2 are two spherical harmonic tensors share the same centre.
     
    
Encode two spherical harmonic tensors that have distinct centers:  

     EquivalentFeatures.W3jProduct(V1,V2,V3) 
	 
This function converts three spherical harmonic tensors into a rotation-invariant vector. These three tensors may correspond to the two spherical harmonic representations of the density or geometry of two objects with distinct centers, along with the tensor obtained from the subtraction between the coordinates of two objects.

     EquivalentFeatures.W3jProductCToR(E1)
	 
This function converts the output of the complex-formatted produced by W3jProduct into a real-valued format.

     EquivalentFeatures.W3jProduct(V1,V2,V3,n1,n2)
	 
This is a restricted version of Wigner 3J. Only shell with degree n1 from V1 is selected and the three shells for the 3j wigner product with degrees s1,s2,s3 satisfying s1 = n1 and s1 < abs(s2-s3)+n2.

     EquivalentFeatures.W3jProductCToR(InvariantV,n1,n2,d2,d3)

This function converts the result of restricted version of Wigner 3J into a real-valued format. The definition of n1 and n2 is the same as W3jProduct. d2,d3 are the degree of V2,V3.  

e.g.

	  
      using  EquivalentFeatures
      using LinearAlgebra



      EquivalentFeatures.Initial();
      S1=[0.28209479177387814 + 0.0im,0.32024057866853123 - 0.20790429007338984im,-0.7888471542040005 + 0.0im,-0.32024057866853123 - 0.20790429007338984im,0.1919939064150656 - 0.4309075639277064im,-1.1561083338265001 + 0.7505603549415122im,1.2590212271988541 - 0.0im,1.1561083338265001 + 0.7505603549415122im,0.1919939064150656 + 0.4309075639277064im]
      S2=[ 0.28209479177387814 + 0.0im,-0.34356259618321017 + 0.42005430129421506im,0.48849798065320743 + 0.0im,0.34356259618321017 + 0.42005430129421506im,-0.1890184428370038 - 0.9340187976228163im,-0.7680649651957328 + 0.9390690252317576im,-0.14756946227815848 + 0.0im,0.7680649651957328 + 0.9390690252317576im,-0.1890184428370038 + 0.9340187976228163im]
      D3 = [0.06750336452980638,0.6140480570957588,0.6147106620083631] 
	  # S1,S2 are two spherical harmonic tensors representing two objects with distinct centers. D3 represents the arrow between two centers.
	  
      S1_= [0.28209479177387814 + 0.0im,-0.4591493984823001 - 0.01920376046845461im,-0.7010437701912248 + 0.0im,0.4591493984823001 - 0.01920376046845461im,0.6810233009616511 + 0.05706694248221995im,1.4730872169146216 + 0.061611349500355135im,0.7405487972724079 - 0.0im,-1.4730872169146216 + 0.061611349500355135im,0.6810233009616511 - 0.05706694248221995im]
      S2_= [0.28209479177387814 + 0.0im,0.5926183815697879 - 0.06877358557907325im,0.3402048201752961 + 0.0im,-0.5926183815697879 - 0.06877358557907325im,1.1211813340313133 - 0.2637795299760591im,0.9226667866122251 - 0.10707582684485425im,-0.6346265484963751 + 0.0im,-0.9226667866122251 - 0.10707582684485425im,1.1211813340313133 + 0.2637795299760591im]
      D3_= [-0.38631329442688933,-0.254642071890792,0.7385122696373354]
      # S1_, S2_, and D_3 are the rotated versions of S1, S2, and D3, respectively.
      
      
      R3 =  norm(D3,2);
      D3 = D3/R3;
      D3 = vcat(1, D3);
      D3 = EquivalentFeatures.CtoS_Encode(D3,2);
      D3 = EquivalentFeatures.IncreaseDegree(D3,1);
      D3[2:4] = D3[2:4]*R3; 
      D3[5:9] = D3[5:9]*R3*R3;
      R3_ =  norm(D3_,2);
      D3_ = D3_/R3_;
      D3_ = vcat(1, D3_);
      D3_ = EquivalentFeatures.CtoS_Encode(D3_,2);
      D3_ = EquivalentFeatures.IncreaseDegree(D3_,1);
      D3_[2:4] = D3_[2:4]*R3_; 
      D3_[5:9] = D3_[5:9]*R3_*R3_;

      E  = EquivalentFeatures.W3jProduct(D3,S1,S2);
      E_ = EquivalentFeatures.W3jProduct(D3_,S1_,S2_); 
      E  = EquivalentFeatures.W3jProductCToR(E);
      E_ = EquivalentFeatures.W3jProductCToR(E_);
      Loss = sum(abs.(E - E_)) 
      E  = EquivalentFeatures.W3jProduct(D3,S1,S2,2,2)
	    # The restriction is that in D3, only the shell with quantum number *l=1* (for *n1=l+1*,n1=2) is used in the Wigner 3J product.
      E_ = EquivalentFeatures.W3jProduct(D3_,S1_,S2_,2,2)
      E  = EquivalentFeatures.W3jProductCToR(E,2,2,3,3)
      E_ = EquivalentFeatures.W3jProductCToR(E_,2,2,3,3)
      Loss = sum(abs.(E - E_)) 




Encoding in a rotation-equivariant manner:

     EquivalentFeatures.W3jProductCompact(V1,R2,R3,n1);
     
This function converts the spherical harmonic tensor V1 into a rotation-invariant vector using the reference spherical harmonic tensors R2 and R3, where n1 denotes the maximum degree of the shell in V1. The resulting invariance vector has the same length as V1. 
   
     EquivalentFeatures.DecodeMatrixCompact(R2,R3,n1); 

This function will return a matrix that can be used to convert the invariant coding calculated by the W3jProductCompact function into one of the original spherical harmonic tensors (V1).  Because the length of V1 is equal to invariant coding, this matrix is therefore a square matrix. And the R2, R3 are two reference spherical harmonic tensors.

e.g.
      using LinearAlgebra
      using EquivalentFeatures


      C1 =  [0.0043477849927746155,0.0,0.9999905483381614];
      R1 =  norm(C1,2);
      C1 =  C1/R1; # Cartesian coordinate normalization.
      V1 =  vcat(1, C1);
      
      C2 =  [0.772027518982468,0.33454525822573616,0.5404192632877276];
      R2 =  norm(C2,2);
      C2 =  C2/R2;# Cartesian coordinate normalization.
      V2 =  vcat(1, C2);

      C3 =  [0.654676578,0.7654334543987,0.765432216890];
      R3 =  norm(C3,2);
      C3 =  C3/R3;# Cartesian coordinate normalization.
      V3 =  vcat(1, C3);
      

      
      V1 = EquivalentFeatures.CtoS_Encode(V1,2);
      V2 = EquivalentFeatures.CtoS_Encode(V2,2);
      V3 = EquivalentFeatures.CtoS_Encode(V3,2);
      
      V1 = EquivalentFeatures.IncreaseDegree(V1,1);
      V2 = EquivalentFeatures.IncreaseDegree(V2,1);
      V3 = EquivalentFeatures.IncreaseDegree(V3,1);  
      # V1 is the spherical harmonic tensor to be encoded. V2,V3 are two reference spherical harmonic tensors.   
	  
      W3 = EquivalentFeatures.W3jProductCompact(V1,V2,V3,3)
	  # W3 is the invariant encoding derived from the W3j product of these three tensors.
      M  = EquivalentFeatures.DecodeMatrixCompact(V2,V3,3) 
	  # M is the matrix that converts the invariant coding (W3) into the original spherical harmonic tensor (V1),the length of V1 is equal to the length of W3, therefore M is a square matrix. 
      norm(M*W3-V1,2);


  
 
References spherical harmonic tensors extraction:
 

      EquivalentFeatures.ReferencesExtract(V_input)
      
This function will yield two reference vectors from the diople and quadruple shells of V_input. And these two reference vectors are ​converted into spherical harmonic tensors via *CtoS_Encode* and *IncreaseDegree* functions. 

    
      EquivalentFeatures.SelfProductMatrix( V_input, int n, int n2);  
      
      
For the Clebsch–Gordan (CG) product of any two shells in the V_input spherical harmonic tensor, only top n2  with smaller degree are selected. The spherial harmonic tensors with degree n are returned.

ReferencesExtract,SelfProductMatrix can be used to generate the reference vectors for W3jProductCompact and DecodeMatrixCompact functions. 

e.g 

     V_input = [ 0.8462843753216345 + 0.0im, 0.49442005281239043 - 0.38003620975475416im,1.1266402071841877 + 0.0im,-0.49442005281239043 - 0.38003620975475416im,0.12624905371373815 - 0.586665146512387im,0.712812513540789 - 0.5922990542960181im,0.8306662145480092 + 0.0im,-0.712812513540789 - 0.5922990542960181im,0.12624905371373815 + 0.586665146512387im ]
     V_output = [ 0.5641895835477564 + 0.0im,0.046180362258809796-0.5151788499369871im,0.563454914288809,-0.046180362258809796-0.5151788499369871im,-0.35831018662719905 -0.036901839488305416im,0.09441424373213234-0.6546737533363572im,0.008536735434823653+0.0im,-0.09441424373213234-0.6546737533363572im,-0.35831018662719905+0.036901839488305416im]
      using  EquivalentFeatures
      
      MM =  EquivalentFeatures.ReferencesExtract(V_input)
      v1 = vcat(1,MM[:,1]);
      v2 = vcat(1,MM[:,2]);# MM[:,1] and MM[:,2] have been normalize.
      v2 = EquivalentFeatures.RStoCS_Encode(v2,2);
      v1 = EquivalentFeatures.RStoCS_Encode(v1,2);
      v1 = EquivalentFeatures.IncreaseDegree(v1,1);
      v2 = EquivalentFeatures.IncreaseDegree(v2,1);
	  # Extract two reference spherical harmonic tensors from the V_input spherical harmonic tensor. 
      V_Encode = EquivalentFeatures.W3jProductCompact(V_output,v1,v2,3);
      M = EquivalentFeatures.DecodeMatrixCompact(v1,v2,3);
      norm(M*V_Encode - V_output,2)
      
      
      MM=EquivalentFeatures.SelfProductMatrix(V_input,2,2)
      ReciprocalRadii   = (MM[1,2]*MM[1,2]-MM[1,1]*MM[1,3]*2)/(0.25*(3/pi));
      ReciprocalRadii   = (ReciprocalRadii)^(-(1/2));
      v1 = vcat(1,MM[1,:].*ReciprocalRadii);
      ReciprocalRadii   = (MM[2,2]*MM[2,2]-MM[2,1]*MM[2,3]*2)/(0.25*(3/pi));
      ReciprocalRadii   = (ReciprocalRadii)^(-(1/2));
      v2 = vcat(1,MM[2,:].*ReciprocalRadii);   # normalize MM[1,:] and MM[2,:]. 
      v2 = EquivalentFeatures.RStoCS_Encode(v2,2);
      v1 = EquivalentFeatures.RStoCS_Encode(v1,2);
      v1 = EquivalentFeatures.IncreaseDegree(v1,1);
      v2 = EquivalentFeatures.IncreaseDegree(v2,1);
      V_Encode = EquivalentFeatures.W3jProductCompact(V_output,v1,v2,3);
      M = EquivalentFeatures.DecodeMatrixCompact(v1,v2,3)
      norm(M*V_Encode - V_output,2);



Derivatives of spherical harmonics: 

       EquivalentFeatures.DerivativeSH(Y)

This function will return the derivative of a Spherical Harmonic tensor with respect to $\theta$ and $\phi$.

      EquivalentFeatures.DerivativeSH_XYZ(Y,DR,ReciprocalRadii,ReciprocalF)
      
This function returns the derivative of a spherical harmonic tensor with respect to Cartesian coordinates. DR is a vector whose elements are fractions of the derivative of radius function and radius function((df(r)/dr)/f(r)). ReciprocalRadii is the reciprocal of the radius in Cartesian coordinates (1/r), ReciprocalF is the reciprocal of the radius function. 

e.g.
      using  EquivalentFeatures
      
      x=0.31
      y=0.14
      z=0.41
      r = (x^2+y^2+z^2)^0.5;
      z  = (z/r);
      x  = (x/r);
      y  = (y/r);
      V3 =[1,x,y,z]
      S = EquivalentFeatures.CtoS_Encode(V3,2);
      S = EquivalentFeatures.IncreaseDegree(S,1);
      S = EquivalentFeatures.IncreaseDegree(S,1);
      f0 = 1;
      f1 = r;
      f2 = r*r;
      f3 = r*r*r; 
      S[1]     = S[1]*f0;
      S[2:4]   = S[2:4]*f1; 
      S[5:9]   = S[5:9]*f2; 
      S[10:16] = S[10:16]*f3; 
      EquivalentFeatures.DerivativeSH(S);
      D = zeros(16);
      D[2:4] .= 1.0;  #f1'
      D[5:9] .= 2*r;  #f2'
      D[10:16] .= 3*r*r; #f3'
      ReciprocalRadii = 1/r; 
      ReciprocalF     = 1/f1;
      EquivalentFeatures.DerivativeSH(Y,DR,ReciprocalRadii,ReciprocalF)
The DerivativeSH function returns a matrix with two columns. One column indicates the θ and the other indicates the φ.  DerivativeSH_XYZ returns a three-column matrix indicating X, Y, and Z, respectively. 

       EquivalentFeatures.DerivativeWignerD(Y)

This function will return the derivative of a bulk level Spherical Harmonic tensors with respect to $\theta$ and $\phi$ using wigner D matrix.

e.g.

      
      using   EquivalentFeatures
      using   WignerD 
      using   LinearAlgebra


      EquivalentFeatures.Initial();
      C1 =  [0.0043477849927746155,0.0,0.9999905483381614]; 
      C2 =  [0.772027518982468,0.33454525822573616,0.5404192632877276];
      C3 =  [0.0043477849927746155,0.9999905483381614,0.11];


       R  =  norm(C1,2);
       C1_ = C1/R; 
       V1 =  vcat(1, C1_);
       S1 = EquivalentFeatures.CtoS_Encode(V1,2);
       S1 = EquivalentFeatures.IncreaseDegree(S1,1);
       S1[2:4] = S1[2:4]*R; 
       S1[5:9] = S1[5:9]*R*R; 
       Idential = S1; 

       R  =  norm(C2,2);
       C2_ = C2/R; 
       V1 =  vcat(1, C2_);
       S1 = EquivalentFeatures.CtoS_Encode(V1,2);
       S1 = EquivalentFeatures.IncreaseDegree(S1,1);
       S1[2:4] = S1[2:4]*R; 
       S1[5:9] = S1[5:9]*R*R; 
       Idential = Idential+S1; 
	   
       R  =  norm(C3,2);
       C3_ = C3/R; 
       V1 =  vcat(1, C3_);
       S1 = EquivalentFeatures.CtoS_Encode(V1,2);
       S1 = EquivalentFeatures.IncreaseDegree(S1,1);
       S1[2:4] = S1[2:4]*R; 
       S1[5:9] = S1[5:9]*R*R; 
       Idential = Idential+S1; 

       Deriva = EquivalentFeatures.DerivativeWignerD(Idential);

       sum(abs.(Deriva[2:4,1] - ((wignerD(1,0,0.0001,0) - wignerD(1,0,0,0))/0.0001)*Idential[2:4]))
       sum(abs.(Deriva[5:9,1] - ((wignerD(2,0,0.0001,0) - wignerD(2,0,0,0))/0.0001)*Idential[5:9]))

       sum(abs.(Deriva[2:4,2] - ((wignerD(1,0.0001,0,0) - wignerD(1,0,0,0))/0.0001)*Idential[2:4])) 
       sum(abs.(Deriva[5:9,2] - ((wignerD(2,0.00001,0,0) - wignerD(2,0,0,0))/0.00001)*Idential[5:9]))
