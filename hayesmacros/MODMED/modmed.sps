/* MODMED version 3.1, posted online on Jan 21, 2011 */.
/* Written by Andrew F. Hayes */.
/* The Ohio State University */.
/* hayes.338@osu.edu */.



/* Note: The syntax structure differs slightly from the structure described in Preacher, Rucker, and Hayes (2007) */.
/* This version produces bootstrap results only when used in conjunction with the dvmodv or mmodv commands */.


define MODMED (vars = !charend ('/')/dv = !charend('/')/med = !charend('/')/dvmodel = !charend('/')/mmodel = !charend('/')
   /dvmodv = !charend ('/') !default (9999)/mmodv = !charend ('/') !default (9999)/covmat = !charend('/') !default(0)
   /varord = !charend('/') !default(2)/boot = !charend('/') !default (1)/jn = !charend('/') !default(0)).
preserve.
set seed = random.
set length = none.
set mxloop = 10000000.
matrix.
do if (!boot > 999).
  compute btn = trunc(!boot/1000)*1000.
  else.
  compute btn = 1.
end if.
compute bootcool = 0.
compute jnon = !jn.
compute p0=-.322232431088.
compute p1 = -1.
compute p2 = -.342242088547.
compute p3 = -.0204231210245.
compute p4 = -.0000453642210148.
compute q0 = .0993484626060.
compute q1 = .588581570495.
compute q2 = .531103462366.
compute q3 = .103537752850.
compute q4 = .0038560700634.

/* This section brings in the data into different matrices */.


get y/file = */variables = !dv/names = ynm/missing = omit.
get med/file = */variables = !med/names = mnm/missing = omit.
get ymodel/file = */variables = !dvmodel/names = ymnms/missing = omit.
get mmodel/file = */variables = !mmodel/names mmnms/missing = omit.
get vars/file = */variables = !vars/names = vname/missing = omit.
compute n = nrow(vars).
compute nvarch = make(1,ncol(vars),0).

/* extract the dv */.
loop i = 1 to ncol(vname).
do if (vname(:,i)=ynm).
  compute y = vars(:,i).
  compute nvarch(1,i)=1.
end if.
end loop.

/* extract the mediator */.
loop i = 1 to ncol(vname).
do if (vname(:,i)=mnm).
  compute med = vars(:,i).
  compute nvarch(1,i)=1.
end if.
end loop.

/* extract the ymodel variables */.
compute ymodel = make(n,ncol(ymnms),0).
loop j = 1 to ncol(ymnms).
loop i = 1 to ncol(vname).
do if (vname(:,i)=ymnms(:,j)).
    compute ymodel(:,j)=vars(:,i).
    compute nvarch(1,i)=1.
end if.
end loop.
end loop.

/* extract the mmodel variables */.
compute mmodel = make(n,ncol(mmnms),0).
loop j = 1 to ncol(mmnms).
loop i = 1 to ncol(vname).
do if (vname(:,i)=mmnms(:,j)).
    compute mmodel(:,j)=vars(:,i).
    compute nvarch(1,i)=1.
end if.
end loop.
end loop.

/* extract the covariates */.
compute vsel = rsum(nvarch).
compute numcovs = ncol(vars)-vsel.
do if (numcovs > 0).
compute covvar = make(n,numcovs,0).
compute covnms = {"x"}.
compute j = 1.
loop i = 1 to ncol(vname).
  do if (nvarch(1,i)) = 0.
  compute covvar(:,j) = vars(:,i).
  compute nvarch(1,i)=1.
  compute j=j+1.
  compute covnms = {covnms, vname(:,i)}.
  end if.
end loop.
compute covnms = covnms(1,2:ncol(covnms)).
end if.

compute bootr = make(btn+1+n,1,0).
compute con = make(n,1,1).
compute varo = !varord.
do if (varo <> 2).
   compute varo = 1.
end if.
print/title = "=======================================================".
print/title = "Preacher, Rucker, & Hayes Moderated Mediation Analysis".
print/title = "=======================================================".
/* Here we do some syntax error checks */.

compute check = ncol(y).
do if (check > 1).
  print/title = "Error: Too many outcome variables in your syntax".
end if.
compute check = ncol(med).
do if (check > 1).
  print/title = "Error: Too many mediator variables in your syntax".
end if.
compute check = ncol(mmodel).
do if (check > 2).
  print/title = "Error: Too many variables in your mediator model".
end if.
compute check = ncol(ymodel).
do if (check > 2).
  print/title = "Error: Too many variables in your y model".
end if.

/* Here we figure out which model is being specified and whether moderators are dichotomous */.
compute mmvar = ncol(mmodel).
compute ymvar = ncol(ymodel).
compute mmdich = 0.
compute ymdich = 0.
do if (mmvar = 2).
compute mmdich = (ncol(design(mmodel(:,2))) = 2).
end if.
do if (ymvar = 2).
compute ymdich = (ncol(design(ymodel(:,2))) = 2).
end if.
compute model = 0.
compute test3 = csum(ymodel(:,1)=med).
compute test = csum(mmodel(:,mmvar)=ymodel(:,ymvar)).
compute test2 = csum(mmodel(:,1)=ymodel(:,1)).
do if (mmvar = 2 AND ymvar = 2 and test <> n and test2 <> n).
   compute model = 4.
   compute interact =mmodel(:,1)&*mmodel(:,2).
   do if (numcovs > 0).
   compute dmdl = {con,mmodel,interact, covvar}.
   compute ivnms = t({"Constant", mmnms, "Inter1", covnms}).
   compute inter2 = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con, mmodel, interact,ymodel,inter2, covvar}.
   compute ivnms2 = t({"Constant", mmnms, "Inter1", ymnms, "Inter2", covnms}).
   end if.
   do if (numcovs = 0).
   compute dmdl = {con,mmodel,interact}.
   compute ivnms = t({"Constant", mmnms, "Inter1"}).
   compute inter2 = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con, mmodel, interact,ymodel,inter2}.
   compute ivnms2 = t({"Constant", mmnms, "Inter1", ymnms, "Inter2"}).
   end if.
   compute prod = {mmnms(1,1), "*", mmnms(1,2); ymnms(1,1), "*", ymnms(1,2)}.
   compute cprod = {"Inter1:  "; "Inter2:  "}.
   compute mv = {mmnms(1,2), ymnms(1,2)}.
   compute mat = csum(mmodel(:,2))/n. 
   compute ones = make(n,1,1).
   compute xyc = mmodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (mmdich = 0).
   compute matm = {mat-sd; mat; mat+sd}.
   end if.
   do if (mmdich = 1).
   compute matm = {cmin(mmodel(:,2)); cmax(mmodel(:,2))}.
   end if. 
   compute mat = csum(ymodel(:,2))/n.
   compute xyc = ymodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (ymdich = 0).
     compute maty = {mat-sd; mat; mat+sd}.
   end if.
   do if (ymdich = 1).
     compute maty = {cmin(ymodel(:,2)); cmax(ymodel(:,2))}.
   end if.
   compute spe = csum(mmodel(:,2))/n.
   compute spe2 = csum(ymodel(:,2))/n. 
   do if (!mmodv = 9999 and !dvmodv = 9999).
       compute mmm = spe.
       compute yyy = spe2.
   end if.
   do if (!mmodv = 9999 and !dvmodv <> 9999).
      compute maty = !dvmodv.
      compute mmm = spe.
      compute yyy = !dvmodv.
   end if.
   do if (!mmodv <> 9999 and !dvmodv = 9999).
      compute matm = !mmodv.
      compute mmm = !mmodv.
      compute yyy = spe2.
   end if.
   do if (!mmodv <> 9999 and !dvmodv <> 9999).
      compute matm = !mmodv.
      compute maty = !dvmodv.
      compute mmm = !mmodv.
      compute yyy = !dvmodv.
   end if.
   compute mat = make((nrow(matm)*nrow(maty)),2,0).
   compute k = 1.
   loop i = 1 to nrow(matm).
     loop j = 1 to nrow(maty).
     compute mat(k,:) = {matm(i,1), maty(j,1)}.
     compute k=k+1.
     end loop.
   end loop.
   compute rmat1 = nrow(mat).
   compute inmnm = mmnms(1,2).
   compute inmnm2 = ymnms(1,2).
   compute cmmin = mmin(mmodel(:,2)).
   compute cmmax = mmax(mmodel(:,2)).
   compute cmmin2 = mmin(ymodel(:,2)).
   compute cmmax2 = mmax(ymodel(:,2)).
else if (mmvar = 2 AND ymvar = 2 AND test = n and test2 <> n).
   compute model = 5.
   compute interact = mmodel(:,1)&*mmodel(:,2).
   do if (numcovs > 0).
   compute dmdl = {con,mmodel,interact, covvar}.
   compute ivnms = t({"Constant", mmnms, "Inter1", covnms}).
   compute interac2 = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con,mmodel,interact,med,interac2, covvar}.
   compute ivnms2 = t({"Constant", mmnms, "Inter1", mnm, "Inter2", covnms}).
   end if.
   do if (numcovs = 0).
   compute dmdl = {con,mmodel,interact}.
   compute ivnms = t({"Constant", mmnms, "Inter1"}).
   compute interac2 = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con,mmodel,interact,med,interac2}.
   compute ivnms2 = t({"Constant", mmnms, "Inter1", mnm, "Inter2"}).
   end if.
   compute prod = {mmnms(1,1), "*", mmnms(1,2); ymnms(1,1), "*", ymnms(1,2)}.
   compute cprod = {"Inter1:  "; "Inter2:  "}.
   compute mv = {mmnms(1,2), " "}.
   compute mat = csum(mmodel(:,2))/n. 
   compute ones = make(n,1,1).
   compute xyc = mmodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (mmdich = 0).
   compute mat = {mat-sd; mat; mat+sd}.
   end if.
   do if (mmdich = 1).
   compute mat = {cmin(ymodel(:,2));cmax(ymodel(:,2))}.
   end if.
   do if (!mmodv <> 9999).
        compute mat = !mmodv.
   end if.
   do if (!dvmodv <> 9999).
        compute mat = !dvmodv.
   end if.
   compute rmat1 = nrow(mat).
   compute cmmin = mmin(mmodel(:,2)).
   compute cmmax = mmax(mmodel(:,2)).   
   compute inmnm = mmnms(1,2).
else if (mmvar = 2 and ymvar = 1 and test <> n and test2 <> n).
   compute model = 2.
   compute interact = mmodel(:,1)&*mmodel(:,2).
   do if (numcovs > 0).
   compute dmdl = {con,mmodel, interact, covvar}.
   compute ivnms = t({"Constant", mmnms, "Inter1", covnms}).
   compute dymdl = {con, ymodel, mmodel, interact, covvar}.
   compute ivnms2 = t({"Constant", ymnms, mmnms, "Inter1", covnms}).
   end if.
   do if (numcovs = 0).
   compute dmdl = {con,mmodel, interact}.
   compute ivnms = t({"Constant", mmnms, "Inter1"}).
   compute dymdl = {con, ymodel, mmodel, interact}.
   compute ivnms2 = t({"Constant", ymnms, mmnms, "Inter1"}).
   end if.
   compute prod = {mmnms(1,1), "*", mmnms(1,2)}.
   compute cprod = {"Inter1: "}.
   compute mv = {mmnms(1,2), " "}.
   compute mat = csum(mmodel(:,2))/n. 
   compute ones = make(n,1,1).
   compute xyc = mmodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (mmdich = 0).
   compute mat = {mat-sd; mat; mat+sd}.
   end if.
   do if (mmdich = 1).
   compute mat = {cmin(mmodel(:,2));cmax(mmodel(:,2))}.
   end if.
   do if (!mmodv <> 9999).
     compute mat=!mmodv.
   end if.
   compute rmat1 = nrow(mat).
   compute cmmin = mmin(mmodel(:,2)).
   compute cmmax = mmax(mmodel(:,2)).
   compute inmnm = mmnms(1,2).
else if (mmvar = 1 and ymvar = 2 and test2 <> n and test <> n).
   compute model = 3.
   do if (numcovs > 0).
   compute dmdl = {con,mmodel, covvar}.
   compute ivnms = t({"Constant", mmnms, covnms}).
   compute interact = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con, mmodel, ymodel, interact, covvar}.
   compute ivnms2 = t({"Constant", mmnms, ymnms, "Inter2", covnms}).
   end if.
   do if (numcovs = 0).
   compute dmdl = {con,mmodel}.
   compute ivnms = t({"Constant", mmnms}).
   compute interact = ymodel(:,1)&*ymodel(:,2).
   compute dymdl = {con, mmodel, ymodel, interact}.
   compute ivnms2 = t({"Constant", mmnms, ymnms, "Inter2"}).
   end if.
   compute prod = {ymnms(1,1), "*", ymnms(1,2)}.
   compute cprod = {"Inter2: "}.
   compute mv = {ymnms(1,2), " "}.
   compute mat = csum(ymodel(:,2))/n. 
   compute ones = make(n,1,1).
   compute xyc = ymodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (ymdich = 0).
   compute mat = {mat-sd; mat; mat+sd}.
   end if.
   do if (ymdich = 1).
   compute mat = {cmin(ymodel(:,2)); cmax(ymodel(:,2))}.
   end if.
   do if (!dvmodv <> 9999).
     compute mat=!dvmodv.
   end if.
   compute rmat1 = nrow(mat).
   compute cmmin = mmin(ymodel(:,2)).
   compute cmmax = mmax(ymodel(:,2)).
   compute inmnm = ymnms(1,2).
else if (mmvar = 1 and ymvar = 2 and test2 = n and test <> n).
  compute model = 0.
else if (mmvar = 1 and ymvar = 2 and test = n and test2 <> n).
   compute model = 1.
   do if (numcovs > 0).
   compute dmdl = {con,mmodel, covvar}.
   compute ivnms = t({"Constant", mmnms, covnms}).
   compute interact = ymodel(:,1)&*ymodel(:,2).
   compute ivnms2 = t({"Constant", ymnms, "Inter2", covnms}).
   compute dymdl = {con,ymodel, interact, covvar}.
   end if.
   do if (numcovs = 0).
   compute dmdl = {con,mmodel}.
   compute ivnms = t({"Constant", mmnms}).
   compute interact = ymodel(:,1)&*ymodel(:,2).
   compute ivnms2 = t({"Constant", ymnms, "Inter2"}).
   compute dymdl = {con,ymodel, interact}.
   end if.
   compute prod = {ymnms(1,1), "*", ymnms(1,2)}.
   compute cprod = {"Inter2: "}.
   compute mv={mmnms(1,1), " "}.
   compute mat = csum(ymodel(:,2))/n. 
   compute ones = make(n,1,1).
   compute xyc = ymodel(:,2)-(ones*mat).
   compute sd = sqrt((1/(n-1))*(t(xyc)*xyc)).
   do if (ymdich = 0).
   compute mat = {mat-sd; mat; mat+sd}.
   end if.
   do if (ymdich = 1).
   compute mat = {cmin(ymodel(:,2)); cmax(ymodel(:,2))}.
   end if.
   do if (!dvmodv <> 9999).
     compute mat=!dvmodv.
   end if.
   compute rmat1 = nrow(mat).
   compute cmmin = mmin(ymodel(:,2)).
   compute cmmax = mmax(ymodel(:,2)).
   compute inmnm = ymnms(1,2).
end if.

do if (test3 <> n).
   compute model = 0.
end if.
do if (model = 0).
   print/title = "It appears you have not specified any model correctly with your syntax".
end if.

do if (model <> 4 and ((!mmodv <> 9999) or (!dvmodv <> 9999))).
  compute jnon = 0.
end if.


/* Here we start the computations */.
do if (numcovs > 0).
compute blks = {" "}.
loop bkl = 1 to ncol(covvar).
compute blks = {blks; " "}.
end loop.
compute blks = blks(2:(ncol(covvar)+1),1).
end if.

do if (model > 0).
print model/title = "You specified model number:".
do if (numcovs > 0).
compute vs = {mmnms(1,1), " "; ynm, " "; mnm, " "; (mv); t(covnms), blks}.
print vs/title = "Variables in System:"/rlabels = "     IV:" "     DV:" "Med Var:" "Mod Var:" "   Covs:"/format a8.
end if.
do if (numcovs = 0).
compute vs = {mmnms(1,1), " "; ynm, " "; mnm, " "; (mv)}.
print vs/title = "Variables in System:"/rlabels = "     IV:" "     DV:" "Med Var:" "Mod Var:"/format a8.
end if.
print n/title = "Sample size:".
compute bone = inv(t(dmdl)*dmdl)*t(dmdl)*med).
compute msres = csum((med-(dmdl*bone))&**2)/(n-nrow(bone)).
compute cone=msres*inv(t(dmdl)*dmdl).
compute sebone = sqrt(diag(cone)).
compute te = bone&/sebone.
compute p = 2*(1-tcdf(abs(te), n-nrow(bone))).
compute oput = {bone,sebone, te, p}.
print/title = '----------'.
print oput/title = 'MEDIATOR VARIABLE MODEL'/clabels = "Coeff" "SE" "t" "P>|t|"/rnames = ivnms/format f10.4.
compute btwo = inv(t(dymdl)*dymdl)*t(dymdl)*y).
compute msres = csum((y-(dymdl*btwo))&**2)/(n-nrow(btwo)).
compute ctwo=msres*inv(t(dymdl)*dymdl).
compute sebtwo = sqrt(diag(ctwo)).
compute te = btwo&/sebtwo.
compute p = 2*(1-tcdf(abs(te), n-nrow(btwo))).
compute oput = {btwo,sebtwo, te, p}.
do if (!covmat = 1).
  print cone/title = "Variance-Covariance Matrix of Estimates"/cnames = ivnms/rnames = ivnms/format F8.4. 
end if.
print/title = '----------'.
print oput/title = 'DEPENDENT VARIABLE MODEL'/clabels = "Coeff" "SE" "t" "P>|t|"/rnames = ivnms2/format f10.4.
do if (!covmat = 1).
  print ctwo/title = "Variance-Covariance Matrix of Estimates"/cnames = ivnms2/rnames = ivnms2/format F8.4. 
end if.
print/title = '----------'.
print prod/title = "Interaction Terms:"/rnames = cprod/format = a8.

/* Here we calculate the ends of the regions of significance */.
   compute cs1 = cmmin-1.
   compute cs2 = cmmin-1.
do if (model = 1 and ymdich = 0).
   compute a12 = bone(2,1)*bone(2,1).
   compute b22 = btwo(4,1)*btwo(4,1).
   compute sa12 = cone(2,2).
   compute sb22 = ctwo(4,4).
   compute a1 = bone(2,1).
   compute b2 = btwo(4,1).
   compute b1 = btwo(2,1).
   compute cb1b2 = ctwo(4,2).
   compute sb12 = ctwo(2,2).
   compute b12 = btwo(2,1)*btwo(2,1).
   compute t1 = 2*((-(a12)*b22)+(b22*sa12*3.8416)+(a12*sb22*3.8416)).
   compute t2 = (2*a12*b1*b2)-(2*a12*cb1b2*3.8416)-(2*b1*b2*sa12*3.8416).
   compute t3 = ((-2*a12*b1*b2)+(2*a12*cb1b2*3.8416)+(2*b1*b2*sa12*3.8416))**2.
   compute t4 = 4*((-(a12)*b12)+(b12*sa12*3.8416)+(a12*sb12*3.8416)).
   compute t5 = (-(a12)*b22)+(b22*sa12*3.8416)+(a12*sb22*3.8416).
   do if (varo = 2).
   compute t1 = t1+(2*(sa12*sb22*3.8416)).
   compute t2 = t2-(2*cb1b2*sa12*3.8416).
   compute t3 = ((-2*a12*b1*b2)+(2*a12*cb1b2*3.8416)+(2*b1*b2*sa12*3.8416)+(2*cb1b2*sa12*3.8416))**2.
   compute t4 = t4+(4*(sa12*sb12*3.8416)).
   compute t5 = t5+(sa12*sb22*3.8416).
   end if.
   compute rad = t3-(t4*t5).
   do if (rad >= 0 and t1 <> 0).
   compute cs1 = (1/t1)* (t2-sqrt(rad)).
   compute cs2 = (1/t1)* (t2+sqrt(rad)).
   end if.
end if.
do if (model = 2 and mmdich = 0).
   compute a32 = bone(4,1)*bone(4,1).
   compute b12 = btwo(2,1)*btwo(2,1).
   compute sa32 = cone(4,4).
   compute sb12 = ctwo(2,2).
   compute a1 = bone(2,1).
   compute a3 = bone(4,1).
   compute ca1a3 = cone(4,2).
   compute a12 = bone(2,1)*bone(2,1).
   compute sa12 = cone(2,2).
   compute t1 = 2*((-(a32)*b12)+(b12*sa32*3.8416)+(a32*sb12*3.8416)).
   compute t2 = (2*a1*a3*b12)-(2*b12*ca1a3*3.8416)-(2*a1*a3*sb12*3.8416).
   compute t3 = ((-2*a1*a3*b12)+(2*b12*ca1a3*3.8416)+(2*a1*a3*sb12*3.8416))**2.
   compute t4 = 4*((-(a12)*b12)+(b12*sa12*3.8416)+(a12*sb12*3.8416)).
   compute t5 = (-(a32)*b12)+(b12*sa32*3.8416)+(a32*sb12*3.8416).
   do if (varo = 2).
   compute t1 = t1+(2*(sa32*sb12*3.8416)).
   compute t2 = t2-(2*ca1a3*sb12*3.8416).
   compute t3 = ((-2*a1*a3*b12)+(2*b12*ca1a3*3.8416)+(2*a1*a3*sb12*3.8416)+(2*ca1a3*sb12*3.8416))**2.
   compute t4 = t4+(4*(sa12*sb12*3.8416)).
   compute t5 = t5+(sa32*sb12*3.8416).
   end if.
   compute rad = t3-(t4*t5).
   do if (rad >= 0 and t1 <> 0).
   compute cs1 = (1/t1)* (t2-sqrt(t3-(t4*t5))).
   compute cs2 = (1/t1)* (t2+sqrt(t3-(t4*t5))).
   end if.
end if.
do if (model = 3 and ymdich = 0).
  compute a12 = bone(2,1)*bone(2,1).
  compute b32 = btwo(5,1)*btwo(5,1).
  compute sa12 = cone(2,2).
  compute sb32 = ctwo(5,5).
  compute b1 = btwo(3,1).
  compute b3 = btwo(5,1).
  compute cb1b3 = ctwo(5,3).
  compute b12 = btwo(3,1)*btwo(3,1).
  compute sb12 = ctwo(3,3).
  compute t1 = 2*((-(a12)*b32)+(b32*sa12*3.8416)+(a12*sb32*3.8416)).
  compute t2 = (2*a12*b1*b3)-(2*a12*cb1b3*3.8416)-(2*b1*b3*sa12*3.8416).
  compute t3 = ((-2*a12*b1*b3)+(2*a12*cb1b3*3.8416)+(2*b1*b3*sa12*3.8416))**2.
  compute t4 = 4*((-(a12)*b12)+(b12*sa12*3.8416)+(a12*sb12*3.8416)).
  compute t5 = (-(a12)*b32)+(b32*sa12*3.8416)+(a12*sb32*3.8416).
  do if (varo = 2).
  compute t1 = t1+(2*(sa12*sb32*3.8416)).
  compute t2 = t2-(2*cb1b3*sa12*3.8416).
  compute t3 = ((-2*a12*b1*b3)+(2*a12*cb1b3*3.8416)+(2*b1*b3*sa12*3.8416)+(2*cb1b3*sa12*3.8416))**2.
  compute t4 = t4+(4*(sa12*sb12*3.8416)).
  compute t5 = t5+(sa12*sb32*3.8416).
  end if.
  compute rad = t3-(t4*t5).
   do if (rad >= 0 and t1 <> 0).
   compute cs1 = (1/t1)* (t2-sqrt(t3-(t4*t5))).
   compute cs2 = (1/t1)* (t2+sqrt(t3-(t4*t5))).
   end if.
end if.
do if (model = 4).
   compute cs1 = cmmin-1.
   compute cs2 = cmmin-1.
   compute cs3 = cmmin2-1.
   compute cs4 = cmmin2-1.
   compute w = mmm.
   compute z = yyy.
   compute a12 = bone(2,1)*bone(2,1).
   compute b32 = btwo(7,1)*btwo(7,1).
   compute sa12 = cone(2,2).
   compute sb32 = ctwo(7,7).
   compute a1 = bone(2,1).
   compute a3 = bone(4,1).
   compute ca1a3 = cone(4,2).
   compute a32 = bone(4,1)*bone(4,1).
   compute sa32 = cone(4,4).
   compute b1 = btwo(5,1).
   compute b3 = btwo(7,1).
   compute cb1b3 = ctwo(7,5).
   compute b12 = btwo(5,1)*btwo(5,1).
   compute sb12 = ctwo(5,5).
   compute a = (-(a12)*b32)+(b32*sa12*3.8416)+(a12*sb32*3.8416)-(2*a1*a3*b32*w)+(2*b32*ca1a3*3.8416*w)+(2*a1*a3*sb32*3.8416*w)-(a32*b32*w*w).
   compute a = a+(b32*sa32*3.8416*w*w)+(a32*sb32*3.8416*w*w).
   compute b = (-2*a12*b1*b3)+(2*a12*cb1b3*3.8416)+(2*b1*b3*sa12*3.8416)-(4*a1*a3*b1*b3*w)+(4*b1*b3*ca1a3*3.8416*w)+(4*a1*a3*cb1b3*3.8416*w).
   compute b = b+(-2*a32*b1*b3*w*w)+(2*a32*cb1b3*3.8416*w*w)+(2*b1*b3*sa32*3.8416*w*w).
   compute c = (-(a12)*b12)+(b12*sa12*3.8416)+(a12*sb12*3.8416)-(2*a1*a3*b12*w)+(2*b12*ca1a3*3.8416*w)+(2*a1*a3*sb12*3.8416*w)-(a32*b12*w*w).
   compute c = c+(b12*sa32*3.8416*w*w)+(a32*sb12*3.8416*w*w).
   do if (varo = 2).
   compute a = a + (sa12*sb32*3.8416)+(2*ca1a3*sb32*3.8416*w)+(sa32*sb32*3.8416*w*w).
   compute b = b + (4*ca1a3*cb1b3*3.8416*w)+(2*cb1b3*sa32*3.8416*w*w)+(2*cb1b3*sa12*3.8416).
   compute c = c + (sa12*sb12*3.8416)+(2*ca1a3*sb12*3.8416*w)+(sa32*sb12*3.8416*w*w).
   end if.
   compute rad = (b*b)-(4*a*c).
   do if (rad >= 0 and a <> 0).
   compute cs3 = (1/(2*a))*(-b+sqrt(rad)).
   compute cs4 = (1/(2*a))*(-b-sqrt(rad)).
   end if.
   compute a = (-(a32)*b12)+(b12*sa32*3.8416)+(a32*sb12*3.8416)-(2*a32*b1*b3*z)+(2*a32*cb1b3*3.8416*z)+(2*b1*b3*sa32*3.8416*z)-(a32*b32*z*z).
   compute a = a + (b32*sa32*3.8416*z*z)+(a32*sb32*3.8416*z*z).
   compute b = (-2*a1*a3*b12)+(2*b12*ca1a3*3.8416)+(2*a1*a3*sb12*3.8416)-(4*a1*a3*b1*b3*z)+(4*b1*b3*ca1a3*3.8416*z)+(4*a1*a3*cb1b3*3.8416*z).
   compute b = b + (-2*a1*a3*b32*z*z)+(2*b32*ca1a3*3.8416*z*z)+(2*a1*a3*sb32*3.8416*z*z).
   compute c = (-(a12)*b12)+(b12*sa12*3.8416)+(a12*sb12*3.8416)-(2*a12*b1*b3*z)+(2*a12*cb1b3*3.8416*z)+(2*b1*b3*sa12*3.8416*z)-(a12*b32*z*z).
   compute c = c + (b32*sa12*3.8416*z*z)+(a12*sb32*3.8416*z*z).
   do if (varo = 2).
   compute a = a + (sa32*sb12*3.8416)+(2*cb1b3*sa32*3.8416*z)+(sa32*sb32*3.8416*z*z).
   compute b = b + (2*ca1a3*sb12*3.8416)+(4*ca1a3*cb1b3*3.8416*z)+(2*ca1a3*sb32*3.8416*z*z).
   compute c = c + (sa12*sb12*3.8416)+(2*cb1b3*sa12*3.8416*z)+(sa12*sb32*3.8416*z*z).
   end if.
   compute rad = (b*b)-(4*a*c).
   do if (rad >= 0 and a <> 0).
   compute cs1 = (1/(2*a))*(-b+sqrt(rad)).
   compute cs2 = (1/(2*a))*(-b-sqrt(rad)).
   end if.
end if.

/* Here we calculate the conditional indirect effects and second order variances */.

do if (model <> 4).
  compute rmat1 = nrow(mat).
  loop j = 0 to 20.
   compute mat = {mat; j}.
  end loop.
  compute temp = mat(1:rmat1,1).
  compute mat = cmmin+mat*((cmmax-cmmin)/20).
  compute mat(1:rmat1,1)=temp.
  compute rmat2 = nrow(mat).
  compute cson = {-9999, -9999}.
  do if (btn = 1).
    loop j = 1 to 20.
      do if (mat((rmat1+j),1) < cs1) and (mat((rmat1+j+1),1) > cs1).
        compute mat = {mat(1:(rmat1+j),1);cs1;mat((rmat1+j+1):rmat2,1)}.
        compute cson(1,1) = cs1.
        compute rmat2 = rmat2+1.
      end if.
      do if (mat((rmat1+j),1) < cs2) and (mat((rmat1+j+1),1) > cs2).
        compute mat = {mat(1:(rmat1+j),1);cs2;mat((rmat1+j+1):rmat2,1)}.
        compute cson(1,2) = cs2.
        compute rmat2 = rmat2+1.
      end if.
    end loop.
  end if.
compute rmat2 = nrow(mat).  
end if.
do if (model = 4).
  compute rmat1 = nrow(mat).
  do if (mmdich <> 1 and ymdich <> 1).
  loop i = 1 to 2.
    loop j = 0 to 20.
      do if (i = 1 and mmdich = 0 and ymdich = 0).
        compute temp = {mmm,(cmmin2+j*((cmmax2-cmmin2)/20))}.
      end if.
      do if (i = 1 and mmdich = 1 and ymdich = 0).
        compute temp = {cmin(mmodel(:,2)),(cmmin2+j*((cmmax2-cmmin2)/20))}.
      end if.
      do if (i = 1 and mmdich = 0 and ymdich = 1).
      compute temp = {(cmmin+j*((cmmax-cmmin)/20)),cmin(ymodel(:,2))}.
      end if.
      do if (i = 2 and mmdich = 0 and ymdich = 0).
          compute temp = {(cmmin+j*((cmmax-cmmin)/20)),yyy}.
      end if.
      do if (i = 2 and mmdich = 1 and ymdich = 0).
        compute temp = {cmax(mmodel(:,2)),(cmmin2+j*((cmmax2-cmmin2)/20))}.
     end if.
     do if (i = 2 and mmdich = 0 and ymdich = 1).
      compute temp = {(cmmin+j*((cmmax-cmmin)/20)),cmax(ymodel(:,2))}.
     end if.
      compute mat = {mat; temp}.
    end loop.
  end loop.
  end if.

  compute rmat2 = nrow(mat).
  compute cson = {-9999, -9999}.
  compute cson2 = {-9999, -9999}.  
  compute pjn1 = 0.
  compute pjn2 = 0.
  do if (btn = 1 and mmdich = 0 and ymdich = 0).
    loop j = 1 to 20.
      do if (mat((rmat1+j),2) < cs3) and (mat((rmat1+j+1),2) > cs3).
        compute mat = {mat(1:(rmat1+j),:);mmm, cs3;mat((rmat1+j+1):rmat2,:)}.
        compute cson2(1,1) = cs3.
        compute pjn1 = 1.
        compute rmat2 = rmat2+1.
      end if.
    do if (mat((rmat1+j),2) < cs4) and (mat((rmat1+j+1),2) > cs4).
       compute mat = {mat(1:(rmat1+j),:);mmm, cs4;mat((rmat1+j+1):rmat2,:)}.
       compute cson2(1,2) = cs4.
       compute pjn2 = 1.
       compute rmat2 = rmat2+1.
    end if.
  end loop.
  do if (rsum(cson2) <> -19998).
  compute temp={-9999}.
  loop j = 1 to 2.
  do if (cson2(1,j) <> -9999).
  compute temp = {temp,cson2(1,j)}.
  end if.
  end loop.
  compute cson2 = temp(1,2:ncol(cson2)).
  end if.
  compute rmat2 = nrow(mat).
  compute cson = {-9999, -9999}.
  loop j = 1 to 20.
    do if (mat((rmat1+20+j),1) < cs1) and (mat((rmat1+j+20+1),1) > cs1).
      compute mat = {mat(1:(rmat1+20+j),:);cs1, yyy;mat((rmat1+j+20+1):rmat2,:)}.
      compute cson(1,2) = cs1.
      compute rmat2 = rmat2+1.
    end if.
    do if (mat((rmat1+20+j),1) < cs2) and (mat((rmat1+j+20+1),1) > cs2).
      compute mat = {mat(1:(rmat1+20+j),:);cs2, yyy;mat((rmat1+j+20+1):rmat2,:)}.
      compute cson(1,2) = cs2.
      compute rmat2 = rmat2+1.
    end if.
  end loop.
  end if.
end if.
do if (rsum(cson) <> -19998).
compute temp={-9999}.
loop j = 1 to 2.
do if (cson(1,j) <> -9999).
compute temp = {temp,cson(1,j)}.
end if.
end loop.
compute cson = temp(1,2:ncol(cson)).
end if.

compute indsum = make(rmat2,1,0).
compute indssq = indsum.

/* Here we do the bootstrapping */.

compute nnn = 0.
do if (btn = 1).
  compute nnn = -n-1.
end if.
loop j = 1 to (btn+1+n+nnn).
  do if (j > 1).
    do if (j < (n+2)).
      do if (j = 2).
         compute bdmdl = dmdl(2:n,:).
         compute bdymdl = dymdl(2:n,:).
         compute bmed = med(2:n,:).
         compute booty = y(2:n,:).
      else if (j = (n+1)).
         compute bdmdl = dmdl(1:(n-1),:).
         compute bdymdl = dymdl(1:(n-1),:).
         compute bmed = med(1:(n-1),:).
         compute booty = y(1:(n-1),:).
      else.
         compute bdmdl = {dmdl(1:(j-2),:);dmdl((j:n),:)}.
         compute bdymdl = {dymdl(1:(j-2),:);dymdl((j:n),:)}.
         compute bmed = {med(1:(j-2),:);med((j:n),:)}.
         compute booty = {y(1:(j-2),:);y((j:n),:)}.
      end if.
    end if.
    do if (j > (n+1)).
      compute v=trunc(uniform(n,1)*n)+1.
      compute bdmdl = dmdl(v,:).
      compute bdymdl = dymdl(v,:).
      compute bmed = med(v,:).
      compute booty = y(v,:).
    end if.
    compute bbone = inv(t(bdmdl)*bdmdl)*t(bdmdl)*bmed).
    compute bbtwo = inv(t(bdymdl)*bdymdl)*t(bdymdl)*booty).
    do if (model = 1).
      compute bind = bbone(2,1)*(bbtwo(2,1)+bbtwo(4,1)*mat(:,1)).
    else if (model = 2).
      compute bind = (bbone(2,1)+bbone(4,1)*mat(:,1))*bbtwo(2,1).
    else if (model = 3).
      compute bind = bbone(2,1)*(bbtwo(3,1)+bbtwo(5,1)*mat(:,1)).
    else if (model = 4).
      compute bind = (bbone(2,1)+bbone(4,1)*mat(:,1))&*(bbtwo(5,1)+bbtwo(7,1)*mat(:,2)).
    else if (model = 5).
      compute bind = (bbone(2,1)+bbone(4,1)*mat(:,1))&*(bbtwo(5,1)+bbtwo(6,1)*mat(:,1)).
    end if.
    compute bootr(j,1)=bind(1,1).
    do if (j > (n+1)).
    compute indsum = indsum+bind.
    compute indssq = indssq + (bind&**2).
    end if.
  end if.
  do if (j = 1).
    do if (model = 1).
      compute spe2 = mat.
      compute ind = bone(2,1)*(btwo(2,1)+btwo(4,1)*spe2).
      do if (varo = 2).
        compute var=(((btwo(2,1)+btwo(4,1)*spe2)&**2)*cone(2,2))+(((bone(2,1)**2)+cone(2,2))*(ctwo(2,2)+(2*ctwo(4,2)*spe2)+(ctwo(4,4)*(spe2&*spe2)))).
      else if (varo = 1).
        compute var=(((btwo(2,1)+btwo(4,1)*spe2)&**2)*cone(2,2))+((bone(2,1)**2)*(ctwo(2,2)+(2*ctwo(4,2)*spe2)+(ctwo(4,4)*(spe2&*spe2)))).
      end if.
      compute seind = sqrt(var).
      compute disp = spe2.
    else if (model = 2).
      compute spe = mat.
      compute ind = (bone(2,1)+bone(4,1)*spe)*btwo(2,1).
      do if (varo = 2).
        compute var=(((bone(2,1)+bone(4,1)*spe)&**2)*ctwo(2,2))+(((btwo(2,1)**2)+ctwo(2,2))*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe)))).
      else if (varo = 1).
        compute var=(((bone(2,1)+bone(4,1)*spe)&**2)*ctwo(2,2))+((btwo(2,1)**2)*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe)))).
      end if.
      compute seind = sqrt(var).
      compute disp = spe.
   else if (model = 3).
      compute spe2 = mat.
      compute ind = bone(2,1)*(btwo(3,1)+btwo(5,1)*spe2).
      do if (varo = 2).
        compute var=(((btwo(3,1)+btwo(5,1)*spe2)&**2)*cone(2,2))+(((bone(2,1)**2)+cone(2,2))*(ctwo(3,3)+(2*ctwo(5,3)*spe2)+(ctwo(5,5)*(spe2&*spe2)))).
      else if (varo = 1).
        compute var=(((btwo(3,1)+btwo(5,1)*spe2)&**2)*cone(2,2))+((bone(2,1)**2)*(ctwo(3,3)+(2*ctwo(5,3)*spe2)+(ctwo(5,5)*(spe2&*spe2)))).
      end if.
      compute seind = sqrt(var).
      compute disp = spe2.
  else if (model = 4).
      compute spe = mat(:,1).
      compute spe2 = mat(:,2).
      compute ind = (bone(2,1)+bone(4,1)*spe)&*(btwo(5,1)+btwo(7,1)*spe2).
      do if (varo = 2).
        compute var = ((bone(2,1)+bone(4,1)*spe)&**2)&*(ctwo(5,5)+(2*ctwo(7,5)*spe2)+(ctwo(7,7)*(spe2&*spe2))).
        compute var = var + ((btwo(5,1)+btwo(7,1)*spe2)&**2)&*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe))).
        compute var = var +  ((ctwo(5,5)+(2*ctwo(7,5)*spe2)+(ctwo(7,7)*(spe2&*spe2)))&*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe)))).
      else if (varo = 1).
        compute var = ((bone(2,1)+bone(4,1)*spe)&**2)&*(ctwo(5,5)+(2*ctwo(7,5)*spe2)+(ctwo(7,7)*(spe2&*spe2))).
        compute var = var + ((btwo(5,1)+btwo(7,1)*spe2)&**2)&*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe))).
      end if.
      compute seind = sqrt(var).
      compute disp = spe.
      compute disp2 = spe2.
  else if (model = 5).
      compute spe = mat.
      compute ind = (bone(2,1)+bone(4,1)*spe)&*(btwo(5,1)+btwo(6,1)*spe).
      do if (varo = 2).
        compute var = ((btwo(5,1)+btwo(6,1)*spe)&**2)&*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe))).
        compute var = var + ((bone(2,1)+bone(4,1)*spe)&**2)&*(ctwo(5,5)+(2*ctwo(6,5)*spe)+(ctwo(6,6)*(spe&*spe))).
        compute var = var +  ((cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe)))&*(ctwo(5,5)+(2*ctwo(6,5)*spe)+(ctwo(6,6)*(spe&*spe)))).
      else if (varo = 1).
        compute var = ((btwo(5,1)+btwo(6,1)*spe)&**2)&*(cone(2,2)+(2*cone(4,2)*spe)+(cone(4,4)*(spe&*spe))).
        compute var = var + ((bone(2,1)+bone(4,1)*spe)&**2)&*(ctwo(5,5)+(2*ctwo(6,5)*spe)+(ctwo(6,6)*(spe&*spe))).
      end if.
      compute seind = sqrt(var).
     compute disp = spe.
  end if.
end if.
end loop.
do if (btn > 1).
compute lvout = bootr(2:(n+1),:).
compute tdotm = csum(lvout)/n.
compute tm = (make(n,1,1))*tdotm.
compute topa = csum((((n-1)/n)*(tm-lvout))&**3).
compute bota = 6*sqrt((csum((((n-1)/n)*(tm-lvout))&**2)&**3)).
compute ahat = topa&/bota.
compute bootmn = indsum/btn.
compute bootse = sqrt(((btn*indssq)-(indsum&**2))/(btn*(btn-1))).
compute bootz = bootmn&/bootse.
compute bootp = 2*(1-cdfnorm(abs(bootz))).
compute bootr = bootr((n+2):(btn+n+1),:).
do if (rmat1 = 1).
  compute bootrtmp = bootr.
  compute bootrtmp(GRADE(bootr)) = bootr.
  compute bootr = bootrtmp.
  compute lower95 = bootr(.025*btn,1).
  compute upper95 = bootr(1+.975*btn,1).
  compute #pv = (bootr < ind(1,1)).
  compute #pv = csum(#pv)/btn.
  compute p = #pv.
  do if (#pv > .5).
    compute p = 1-#pv.
  end if.
  compute y=sqrt(-2*ln(p)).
  compute xp=y+((((y*p4+p3)*y+p2)*y+p1)*y+p0)/((((y*q4+q3)*y+q2)*y+q1)*y+q0).
  do if (#pv <= .5).
    compute xp = -xp.
  end if.
  compute zlo = xp + ((xp+-1.96)&/(1-t(ahat)&*(xp+-1.96))).
  compute zup = xp + ((xp+1.96)&/(1-t(ahat)&*(xp+1.96))).
  compute zlo = cdfnorm(zlo).
  compute zup = cdfnorm(zup).
  compute blow = trunc(zlo*(btn+1)).
  compute bhigh = trunc(zup*(btn+1))+1.
  do if (blow < 1).
    compute blow(i,1) = 1.
  end if.
  compute lowbca=bootr(blow,1).
  do if (bhigh > btn).
    compute bhigh = btn.
  end if.
  compute upbca =bootr(bhigh,1).
compute ahat = 0.
  compute zlo = xp + ((xp+-1.96)&/(1-t(ahat)&*(xp+-1.96))).
  compute zup = xp + ((xp+1.96)&/(1-t(ahat)&*(xp+1.96))).
  compute zlo = cdfnorm(zlo).
  compute zup = cdfnorm(zup).
  compute blow = trunc(zlo*(btn+1)).
  compute bhigh = trunc(zup*(btn+1))+1.
  do if (blow < 1).
    compute blow(i,1) = 1.
  end if.
  compute lowbc=bootr(blow,1).
  do if (bhigh > btn).
    compute bhigh = btn.
  end if.
  compute upbc =bootr(bhigh,1).
  compute bt = {lower95, upper95, lowbc, upbc, lowbca, upbca}.
end if.
end if.

print/title = '----------'.
compute z = ind/seind.
compute p = 2*(1-cdfnorm(abs(z))).
do if (model <> 4).
    compute cnm = {inmnm, "Ind Eff", "SE", "Z" , "P>|Z|"}.
    compute indout = {disp, ind,seind,z, p}.
else if (model = 4).
    compute cnm = {inmnm, inmnm2, "Ind Eff", "SE", "Z" , "P>|Z|"}.
    compute indout = {disp, disp2, ind,seind,z, p}.
end if.

compute indout2 = indout(1:rmat1,:).
print indout2/title = "Conditional indirect effect at specific value(s) of the moderator(s)"/cnames = cnm/format F8.4.

do if (rmat1 = 1 and btn > 1).
  print/title = "Bootstrap 95% Confidence Intervals for Conditional Indirect Effect".
  compute bt1=bt(1,1:2).
  print bt1/title = "Percentile"
      /clabels = "Lower" "Upper"/format F8.4.
  compute bt1=bt(1,3:4).
  print bt1/title = "Bias Corrected"
      /clabels = "Lower" "Upper"/format F8.4.
  compute bt1=bt(1,5:6).
  print bt1/title = "Bias Corrected and Accelerated"
      /clabels = "Lower" "Upper"/format F8.4.
  print btn/title = "Number of bootstrap samples:"/format F8.0.
  compute bootcool = 1.
end if.

do if (rmat1 > 2 and model <> 4).
  print/title = "Moderator values listed are the sample mean and +/- 1 SD".
end if.

do if ((rmat1 =3 or rmat1 = 6) and model = 4).
   print/title = "Moderator values for one moderator are the sample mean and +/- 1 SD".
end if.

do if (rmat1 = 9 and model = 4).
   print/title = "Moderator values for both moderators are the sample mean and +/- 1 SD".
end if.

print/title = '----------'.
do if (((mmdich = 0 and ymdich = 0) or (model = 4 and (mmdich = 0 or ymdich = 0))) and jnon = 1).
compute indout2 = indout((rmat1+1):rmat2,:).
do if (rsum(cson) <> -19998 and model < 4).
print cson/title = "Moderator value(s) in range of data defining Johnson-Neyman Significance Region(s)"/format = F8.4.
else if (rsum(cson) = -19998 and model < 4 and (indout2(1,5) > .05)).
print/title = "There are no moderator values in range of data defining Johnson-Neyman Significance Regions".
else if (rsum(cson) = -19998 and model < 4 and (indout2(1,5)) <= .05 and (indout2(21,5) <= .05)).
print/title = "The conditional indirect effect is significant throughout the range of the moderator variable".
end if.
do if (model = 4 and !mmodv <> 9999 and !dvmodv = 9999).
compute indout2 = indout2(1:(21+pjn1+pjn2),:).
end if.
do if (model = 4 and !mmodv = 9999 and !dvmodv <> 9999).
compute indout2 = indout2((21+pjn1+pjn2+1):nrow(indout2),:).
end if.
do if (model = 5) and (jnon = 1).
print/title = "Johnson-Neyman computations not available for model 5".
print/title = '----------'.
end if.
print indout2/title = "Conditional indirect effect at range of values of the moderator(s)"/cnames = cnm/format F8.4.
print/title = '----------'.
end if.

do if (btn = 1).
  do if (varo = 2).
    print/title = "Conditional indirect effect standard errors are second-order estimates.".
  else if (varo = 1).
    print/title = "Conditional indirect effect standard errors are first-order estimates.".
  end if.
end if.

end if.

do if (bootcool = 0 and btn > 1).
   print/title = "Bootstrapping disabled.  Bootstrap results provided only in conjunction with".
   print/title = "the use of the /dvmodv or /mmodv subcommands".
end if.

do if (btn = 1).
   print/title = "Bootstrap confidence intervals are recommended for inference about conditional indirect effects".
end if.

do if (mmdich = 1).
  print/title = "The moderator in the mediator model is dichotomous".
end if.
do if (ymdich = 1).
  print/title = "The moderator in the outcome model is dichotomous".
end if.

end matrix.
restore.
!enddefine.

