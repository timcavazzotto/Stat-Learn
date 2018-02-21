/* Written by Andrew F. Hayes */
/* School of Communication */
/* The Ohio State University */
/* See Preacher, KJ, & Hayes, AF (2004).  SPSS and SAS Procedures for Estimating */
/* Indirect Effects in Simple Mediation Models, in Behavior Research Methods, Instruments */
/* and Computers, 36, 717-731 */
/* Version 2.1 uploaded April 13, 2009 */

%macro sobel(data=,y=,x=,m=,boot=);                                               
                                                                                  
/* READ ACTIVE SAS DATA FILE */                                                   
proc iml;                                                                         
use &data where (&y ^= . & &x ^= . & &m ^= .) ;                                   
read all var {&y &x &m};                                                          
n=nrow(&y);                                                                       
dd=&y||&x||&m;                                                                    
dat=dd;                                                                           
                                                                                  
/* DEFINE NUMBER OF BOOTSTRAP SAMPLES */                                          
if (&boot > 999) then;                                                            
  do;                                                                             
  btn = floor(&boot/1000)*1000;                                                   
  btnp=btn+1;                                                                     
  end;                                                                            
                                                                                  
/* START OF THE LOOP FOR BOOTSTRAPPING */                                         
if (&boot < 1000) then;                                                           
  do;                                                                             
  btn = 1000;                                                                     
  btnp = btn+1;                                                                   
  end;                                                                            
res=j(btnp,1,0);                                                                  
do mm = 1 to btnp;                                                                
                                                                                  
  /* DO THE RESAMPLING OF THE DATA */                                             
  if (mm > 1) then if (&boot > 999) then;                                         
    do;                                                                           
      do nn=1 to n;                                                               
      v = int(uniform(0)*n)+1;                                                    
      dat[nn,1:3]=dd[v,1:3];                                                      
      end;                                                                        
    end;                                                                          
  con=j(n,1,1);                                                                   
                                                                                  
  /* SET UP THE DATA COLUMNS FOR PROCESSING */                                    
  x=dat[,2];                                                                      
  y=dat[,1];                                                                      
  m=dat[,3];                                                                      
  xt=dat-J(n,1)*dat[:,];                                                          
  cv=(xt`*xt)/(n-1);                                                              
  sd=sqrt(diag(cv));                                                              
  r=inv(sd)*cv*inv(sd);  
  r2my = r[3,1]*r[3,1];
  r2xy = r[2,1]*r[2,1]; 
  /* CALCULATE REGRESSION STATISTICS NEEDED TO COMPUTE c-c'  */                   
  /* c-c' is held as variable 'ind' */                                            
  xo=con||x;                                                                      
  bzx=inv(xo`*xo)*xo`*m;                                                          
  bzx=bzx[2,1];                                                                   
  xzo=con||x||m;                                                                  
  byzx2=inv(xzo`*xzo)*xzo`*y;                                                     
  byzx=byzx2[3,1];                                                                
  byxz=byzx2[2,1];                                                                
  ind=bzx*byzx;                                                                   
  res[mm,1]=ind;                                                                  
                                                                                  
  /* GENERATE STATISTICS FOR BARON AND KENNY AND NORMAL SOBEL SECTION OF OUTPUT */
  if (mm = 1) then;                                                               
    do;                                                                           
      sdbzx=(sd[3,3]/sd[2,2])*sqrt((1-(r[3,2]*r[3,2]))/(n-2));                    
      ryi=r[2:3,1];                                                               
      rii=r[2:3,2:3];                                                             
      bi=inv(rii)*ryi;                                                            
      rsq=ryi`*bi;
      r2med = r2my-(rsq-r2xy); 
      sec=sqrt((1-rsq)/(n-3))*sqrt(1/(1-(r[3,2]*r[3,2])));                        
      sdyzx=(sd[1,1]/sd[3,3])*sec;                                                
      sdyxz=(sd[1,1]/sd[2,2])*sec;                                                
      seind=sqrt(((byzx*byzx)*(sdbzx*sdbzx))+((bzx*bzx)*(sdyzx*sdyzx))+           
       ((sdbzx*sdbzx)*(sdyzx*sdyzx)));                                            
      byx = r[2,1]*sd[1,1]/sd[2,2];                                               
      sebyx=(sd[1,1]/sd[2,2])*sqrt((1-(r[2,1]*r[2,1]))/(n-2));                    
      se = sebyx//sdbzx//sdyzx//sdyxz;                                            
      bb= byx//bzx//byzx//byxz;                                                   
      tt = bb/se;                                                                 
      df=j(4,1,n-2);                                                              
      df[3,1]=n-3;                                                                
      df[4,1]=n-3;                                                                
      p=2*(1-probt(abs(tt),df));                                                  
      bw=bb||se||tt||p;                                                           
      tst=ind/seind;                                                              
      pv=2*(1-probnorm(abs(tst)));                                                
      LL95 = ind-1.96*seind;                                                      
      UL95=ind+1.96*seind;                                                        
      op=ind||seind||LL95||UL95||tst||pv;                                         
    end;                                                                          
end;                                                                              
/* END OF BOOTSTRAPPING LOOP */                                                   
                                                                                  
/* COMPUTE MEAN AND STANDARD DEV OF INDIRECT EFFECT ACROSS BOOTSTRAP SAMPLES */   
res=res[2:btnp,1];                                                                
mnbt = sum(res)/btn;                                                              
res=-999//res;                                                                    
                                                                                  
  /* SORT THE BOOTSTRAP ESTIMATES */                                              
  do i=2 to btnp;                                                                 
    ix=res[i,1];                                                                  
    do k =i to 2 by -1;                                                           
      m=k;                                                                        
      if res[k-1,1] > ix then;                                                    
        do;                                                                       
          res[k,1]=res[k-1,1];                                                    
        end;                                                                      
      else;                                                                       
      if res[k-1,1] <= ix then;                                                   
        do;                                                                       
          goto stpit;                                                             
        end;                                                                      
    end;                                                                          
    stpit:                                                                        
    res[m,1]=ix;                                                                  
  end;                                                                            
res=res[2:btnp,1];                                                                
btpt=sum(abs(res) >= abs(op[1,1]))/btn;                                           
if op[1,1] <=0 then;                                                              
  do;                                                                             
    btpo=sum(res <= op[1,1])/btn;                                                 
  end;                                                                            
else;                                                                             
  if op[1,1] >= 0 then;                                                           
    do;                                                                           
      btpo=sum(res >= op[1,1])/btn;                                               
    end;                                                                          
                                                                                  
/* GENERATE BOOTSTRAP CONFIDENCE INTERVAL FOR INDIRECT EFFECT */                  
lower99 = res[.005*btn,1];                                                        
lower95 = res[.025*btn,1];                                                        
upper95 = res[1+.975*btn,1];                                                      
upper99 = res[1+.995*btn,1];                                                      
xt=res-J(btn,1)*res[:,];                                                          
cv=(xt`*xt)/(btn-1);                                                              
se=sqrt(diag(cv));                                                                
bt=mnbt||se||lower95||upper95||lower99||upper99;                                  
                                                                                  
/* GENERATE OUTPUT */                                                             
rn={"b(YX)" "b(MX)" "b(YM.X)" "b(YX.M)"};                                         
cn={"Coeff" "s.e." "t" "Sig(Two)"};                                               
print "DIRECT AND TOTAL EFFECTS";                                                 
print bw [rowname = rn colname = cn format = 9.4];                                
rn={"Sobel"};                                                                     
cn={"Value" "s.e." "LL 95 CI" "UL 95 CI" "z" "Sig(Two)"};                         
print "ESTIMATE AND TEST OF INDIRECT EFFECT";                                     
print op [rowname = rn colname = cn format= 9.4];  
print "FAIRCHILD ET AL. (2009) VARIANCE IN Y ACCOUNTED FOR BY INDIRECT EFFECT";
print r2med [format = 9.4]; 
if (&boot > 999) then;                                                            
  do;                                                                             
    print "BOOTSTRAP RESULTS FOR INDIRECT EFFECT";                                
    rn={"Effect"};                                                                
    cn={"Mean" "s.e." "LL 95 CI" "UL 95 CI" "LL 99 CI" "UL 99 CI"};               
    print bt [rowname = rn colname = cn format = 9.4];                            
    print "NUMBER OF BOOTSTRAP RESAMPLES" btn;                                    
  end;                                                                            
print "SAMPLE SIZE" n;                                                            
quit;                                                                             
%mend sobel;   
