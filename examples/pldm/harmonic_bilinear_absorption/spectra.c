//
//  spectra.c
//  
//
//  Created by FrancescoSegatta on 08/06/17.
//
//
//NOTES
//KI,KII -> to obtain one or the other just FT in different Omega1 plane (positive or negative)
//KIII -> ALONG WHICH AXIS FT??

//IMPROVEMENTS:
//1DPP -> FT along t2 // plot selected resp function
//2DPP -> FT along t2 and look w2 slices then -but should I subtract exponential decay?- (nb -> will I see here the different dependence of KI and KII on t2?)
//ANALYSIS OF 2D MAPS (follow point in time//analyze cuts)

//SELECT CONTRIBUTIONS (SE, GSB, ESA)
//Dipoles as vectors, averaging over dipole orientations -> ongoing
//Account for the pulse shape -> basic implementation: for every trj keep track of which state you have visited during each interaction and muliply by the pulse shape.. but doing so we are not considering the effect in the time..

//TEST -> FT spectron Resp Function VERIFIED


int read_data_file(char*);


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define fs_to_cm_m1 33356.41
#define ps_to_cm_m1 33.35641

int main(){


    int i, j, k,l,length;
    double action;
    double flag1, flag2;
    double **data_in;
    char nome_file_out[400], nome_file_in[400];
    FILE *fp_write, *fp_read, *gnuplotPipe;
    
    
    action=1;
    while (action!=0) {
        
        printf("\nOptions:\n--------------------------------------------\n\n");
        printf("             *** RESP. FUNCT. --> FT ***         \n\n");
        printf("    *) Plot & Save LA data :                        1\n");
        printf("    *) Plot & Save 1DPP data :                      2\n");
        printf("    *) Plot & Save 2DPP data (MOD):                 3\n\n");
        //printf("    *) Change file:                           10\n");
        printf("    *) Exit:                                        0\n");
        scanf("%lf",&action);
        printf("--------------------------------------------\n");
        
        //FT OF LA DATA
        if(action==1){
            
            int n_w;
            int ipn;
            double **signal[2];
            double Dt, Dw, w_min, w_max;
            
            printf("Insert RF file name (expected format time (fs, first line t=0), RE_f, IM_f): \n");
            scanf("%s",nome_file_in);
            //sprintf(nome_file_in,"pldm.1-3");
            
            
            //---------------------------------READ THE FILE---------------------------------//
            length=read_data_file(nome_file_in);
            
            data_in=malloc(3*sizeof(double*));
            for (i=0; i<3; i++) data_in[i]=malloc(length*sizeof(double));
            
            fp_read=fopen(nome_file_in,"r");
            
            for (i=0; i<length; i++) {
                
                //pldm_out --> data[0] is time, data[1] is real, data[2] is imaginary
                fscanf(fp_read,"%lf\t%lf\t%lf\n",&data_in[0][i],&data_in[1][i],&data_in[2][i]);
                
            }
            
            fclose(fp_read);
            
            Dt=data_in[0][1]-data_in[0][0];
            printf("FILE METRICS:\n");
            printf("nt: %.3f\n", (float)length);
            printf("Dt: %.3f\n", Dt);
            
            
            
            //////////////////// POSSIBLE MODIFICATION OF THE DATA ///////////////////
            
            for (i=0; i<length; i++) {
                
                if(i==0) printf("INTRODUCE 'COS' CORRECTION\n\n");
                data_in[1][i]*=cos(M_PI*data_in[0][i]/(2.*data_in[0][length-1]));
                data_in[2][i]*=cos(M_PI*data_in[0][i]/(2.*data_in[0][length-1]));
                
            }
            
            //////////////////////////////////////////////////////////////////////////
            
            printf("!! READ FILES !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            
            
            //---------------------------------PARAMETERS---------------------------------//
            
            //AUTOMATIC w domain
            //freq in 1/ps
            w_max=2.*M_PI/Dt/2.; // */2 because of nyquist
            w_min=2.*M_PI/(Dt*length)*2.; // *2
            
            n_w=length;
            
            //SET w domain
            printf("SET w DOMAIN cm^-1 (min, max, #): \n");
            scanf("%lf\t%lf\t%d", &w_min, &w_max, &n_w);
            w_min=2.*M_PI*w_min/fs_to_cm_m1;
            w_max=2.*M_PI*w_max/fs_to_cm_m1;
            
            //if(w_min<2.*M_PI/(Dt*length)*2.||w_min>=w_max) w_min=2.*M_PI/(Dt*length)*2.;
            //if(w_max>2.*M_PI/Dt/2.||w_max<=w_min) w_max=2.*M_PI/Dt/2.;
            
            Dw=(w_max-w_min)/(n_w-1);
            
            //-------------------------------------------------------------------------------//
            
            
            
            
            //--------------------------------- COMPUTE FT ---------------------------------//
            //TRAPEZOID RULE! for every w integrate over t (0-to-EOF), (RE_f*cos-IM_f*sin,IM_f*cos+RE_f*sin)
            
            signal[0]=malloc(3*sizeof(double*));
            signal[1]=malloc(3*sizeof(double*));
            for (i=0; i<3; i++) {
                signal[0][i]=calloc(n_w,sizeof(double));
                signal[1][i]=calloc(n_w,sizeof(double));
            }
            
            //cicle over the frequencies
            for(k=0; k<n_w; k++){
                
                for (ipn=0; ipn<2; ipn++) {
                    signal[ipn][0][k]=(w_min+Dw*k)*(1.-2.*ipn);
                    
                    //REAL
                    signal[ipn][1][k]=data_in[1][0]-data_in[2][length-1]*sin(signal[ipn][0][k]*Dt*(length-1))+data_in[1][length-1]*cos(signal[ipn][0][k]*Dt*(length-1));
                    signal[ipn][1][k]/=2;
                    //IMAGINARY
                    signal[ipn][2][k]=data_in[2][0]+data_in[1][length-1]*sin(signal[ipn][0][k]*Dt*(length-1))+data_in[2][length-1]*cos(signal[ipn][0][k]*Dt*(length-1));
                    signal[ipn][2][k]/=2;
                    
                    //cicle over the time (exept 0 and last point, already accounted above)
                    for(i=1; i<length-1; i++){
                        
                        //REAL
                        signal[ipn][1][k]+=-data_in[2][i]*sin(signal[ipn][0][k]*Dt*i)+data_in[1][i]*cos(signal[ipn][0][k]*Dt*i);
                        //IMAGINARY
                        signal[ipn][2][k]+=data_in[1][i]*sin(signal[ipn][0][k]*Dt*i)+data_in[2][i]*cos(signal[ipn][0][k]*Dt*i);
                        
                    }
                    
                }
                
                printf("\rCompleted %.2f\%%",((double)k+1)*100./n_w);
            }
            
            printf("!! FT PERFORMED !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            //--------------------------------- SAVE FILE ---------------------------------//
            
            //IF NEEDED SHIFT THE FREQUENCY DOMAIN (E.G. IN SPECTRON THE CARRYING FREQUENCY IS REMOVED IN TIME DOMAIN AND REINTRODUCED AS A SHIFTING IN THE FREQ DOMAIN)
            
            sprintf(nome_file_out,"FT_%s",nome_file_in);
            fp_write=fopen(nome_file_out,"w");
            
            for(k=0; k<n_w; k++) fprintf(fp_write,"%lf\t%lf\t%lf\n",signal[0][0][k]*fs_to_cm_m1/(2.*M_PI),signal[0][1][k]+signal[1][1][k],signal[0][2][k]+signal[1][2][k]);
            
            fclose(fp_write);
            
            printf("!! FILE WRITTEN !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            /*------------------------------ PLOT RSPF AND FT ------------------------------//
            
            gnuplotPipe=popen(" gnuplot -persist", "w");
            
            fprintf(gnuplotPipe, "%s \n", "set term x11 1");
            fprintf(gnuplotPipe, "%s \n", "set size ratio 0.3");
            fprintf(gnuplotPipe, "%s \n", "set multiplot layout 2, 1 ");
            
            fprintf(gnuplotPipe, "%s \n", "set title \"RESPONSE FUNCTION\"");
            fprintf(gnuplotPipe, "%s \n", "set xlabel \"time (ps)\"");
            //fprintf(gnuplotPipe, "plot '%s' u 1:2 w l lc rgb 'red' lw 3 title 'imag', '%s' u 1:(-$3) w l lc rgb 'blue' lw 1 title 'real' \n",nome_file_in,nome_file_in);
            //to visualize also possible modifications...
            fprintf(gnuplotPipe, "plot '-' u 1:2 w l lc rgb 'red' lw 3 title 'imag', '-' u 1:3 w l lc rgb 'blue' lw 1 title 'real' \n");
            for (i=0; i<length; i++) fprintf(gnuplotPipe, "%lf\t%lf\t%lf\n",data_in[0][i],data_in[1][i],data_in[2][i]);
            fprintf(gnuplotPipe,"e\n");
            for (i=0; i<length; i++) fprintf(gnuplotPipe, "%lf\t%lf\t%lf\n",data_in[0][i],data_in[1][i],data_in[2][i]);
            fprintf(gnuplotPipe,"e\n");
            fprintf(gnuplotPipe, "%s \n", "set title \"RESPONSE FUNCTION (first 10%)\"");
            fprintf(gnuplotPipe, "set xrange [0:%lf]\n",data_in[0][length-1]/10);
            fprintf(gnuplotPipe,"replot\n");
            
            
            fprintf(gnuplotPipe, "%s \n", "reset");
            fprintf(gnuplotPipe, "%s \n", "unset multiplot");
            fprintf(gnuplotPipe, "%s \n", "set term x11 2");
            fprintf(gnuplotPipe, "%s \n", "set title \"LINEAR SPECTRUM\"");
            fprintf(gnuplotPipe, "%s \n", "set xlabel \"wavenumbers (cm^{-1})\"");
            fprintf(gnuplotPipe, "plot '%s' u 1:2 w l lc rgb 'black' lw 3 title '' \n",nome_file_out);
            fflush(gnuplotPipe);
            //getchar();
            //getchar();
            
            pclose(gnuplotPipe);
            //-------------------------------------------------------------------------------*/
            
            //FREE SECTION
            for (i=0; i<3; i++) free(data_in[i]);
            for (i=0; i<3; i++) {
                free(signal[0][i]);
                free(signal[1][i]);
            }
            free(data_in);
            free(signal[0]);
            free(signal[1]);
            
        }
    
        
        
        //FT OF 1DPP DATA
        if(action==2){
            
            int nt2, nt3, nw3;
            double Dt2, Dt3, Dw3, w3_min, w3_max;
            double **signal;

            printf("Insert RF file name (expected format t1=0 (fs), t2 (fs), t3 (fs), RE_f, IM_f): \n");
            scanf("%s",nome_file_in);
            //sprintf(nome_file_in,"pldm.1-3");
            
            
            //---------------------------------READ THE FILE---------------------------------//
            length=read_data_file(nome_file_in);
            
            data_in=malloc(5*sizeof(double*));
            for (i=0; i<5; i++) data_in[i]=malloc(length*sizeof(double));
            
            fp_read=fopen(nome_file_in,"r");
            
            nt3=0;
            for (i=0; i<length; i++) {
                fscanf(fp_read,"%lf\t%lf\t%lf\t%lf\t%lf\n",&data_in[0][i],&data_in[1][i],&data_in[2][i],&data_in[3][i],&data_in[4][i]);
                if(data_in[1][i]==data_in[1][0]) nt3++;
            }
            
            fclose(fp_read);
            
            nt2=length/nt3;
            Dt2=data_in[1][nt3]-data_in[1][0];
            Dt3=data_in[2][1]-data_in[2][0];

            printf("FILE METRICS:\n");
            printf("nt2: %d\tnt3: %d\n", nt2, nt3);
            printf("Dt2: %f\tDt3: %f\n\n", Dt2, Dt3);
            
            //////////////////// POSSIBLE MODIFICATION OF THE DATA ///////////////////
            for (i=0; i<length; i++) {
                //data_in[3][i]*=cos(M_PI*data_in[2][i]/(2.*data_in[2][nt3-1]));
                //data_in[4][i]*=-cos(M_PI*data_in[2][i]/(2.*data_in[2][nt3-1]));
                data_in[4][i]*=-1;
            }
            //////////////////////////////////////////////////////////////////////////
            
            printf("!! READ FILES !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            //---------------------------------PARAMETERS---------------------------------//
            
            //AUTOMATIC w domain
            //freq in 1/fs
            w3_max=2.*M_PI/Dt3/2.; // x/2 because of nyquist
            w3_min=2.*M_PI/(Dt3*nt3)*2.; // *2
            
            nw3=nt3;
            
            //SET w domain
            printf("SET w DOMAIN cm^-1 (min, max, #): \n");
            scanf("%lf\t%lf\t%d", &w3_min, &w3_max, &nw3);
            w3_min=2.*M_PI*w3_min/fs_to_cm_m1;
            w3_max=2.*M_PI*w3_max/fs_to_cm_m1;
            
            //if(w_min<2.*M_PI/(Dt*length)*2.||w_min>=w_max) w_min=2.*M_PI/(Dt*length)*2.;
            //if(w_max>2.*M_PI/Dt/2.||w_max<=w_min) w_max=2.*M_PI/Dt/2.;
            
            Dw3=(w3_max-w3_min)/(nw3-1);
            
            //-------------------------------------------------------------------------------//
            
            
            
            
            //--------------------------------- COMPUTE FT ---------------------------------//
            //TRAPEZOID RULE! for every t2 and every w3 integrate over t3 (0-to-EOF), (RE_f*cos-IM_f*sin,IM_f*cos+RE_f*sin)
            
            signal=malloc(4*sizeof(double*));
            for (i=0; i<4; i++) signal[i]=calloc(nw3*nt2,sizeof(double));
            
            //cicle over t2
            for (l=0; l<nt2; l++) {
                
                //cicle over the frequencies
                for(k=0; k<nw3; k++){
                    
                    //t2 time
                    signal[0][l*nw3+k]=data_in[1][l*nt3];
                    //w3 freq
                    signal[1][l*nw3+k]=w3_min+Dw3*k;
                    
                    //REAL
                    signal[2][l*nw3+k]=data_in[3][l*nt3]-data_in[4][(l+1)*nt3-1]*sin(signal[1][l*nw3+k]*Dt3*(nt3-1))+data_in[3][(l+1)*nt3-1]*cos(signal[1][l*nw3+k]*Dt3*(nt3-1));
                    signal[2][l*nw3+k]/=2;
                    //IMAGINARY
                    signal[3][l*nw3+k]=data_in[4][l*nt3]+data_in[3][(l+1)*nt3-1]*sin(signal[1][l*nw3+k]*Dt3*(nt3-1))+data_in[4][(l+1)*nt3-1]*cos(signal[1][l*nw3+k]*Dt3*(nt3-1));
                    signal[3][l*nw3+k]/=2;
                    
                    //cicle over the time t3 (exept 0 and last point, already accounted above)
                    for(i=1; i<nt3-1; i++){
                        
                        //REAL
                        signal[2][l*nw3+k]+=-data_in[4][i+l*nt3]*sin(signal[1][l*nw3+k]*Dt3*i)+data_in[3][i+l*nt3]*cos(signal[1][l*nw3+k]*Dt3*i);
                        //IMAGINARY
                        signal[3][l*nw3+k]+=data_in[3][i+l*nt3]*sin(signal[1][l*nw3+k]*Dt3*i)+data_in[4][i+l*nt3]*cos(signal[1][l*nw3+k]*Dt3*i);
                    
                    }
                    
                    printf("\rCompleted %.2f\%%",((double)(l*nw3+k)+1.)*100./nt2/nw3);
                }
                

            }
            
            printf("!! FT PERFORMED !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            //--------------------------------- SAVE FILE ---------------------------------//
            
            //IF NEEDED SHIFT THE FREQUENCY DOMAIN (E.G. IN SPECTRON THE CARRYING FREQUENCY IS REMOVED IN TIME DOMAIN AND REINTRODUCED AS A SHIFTING IN THE FREQ DOMAIN)
            
            sprintf(nome_file_out,"FT_%s",nome_file_in);
            fp_write=fopen(nome_file_out,"w");
        
            for (l=0; l<nt2; l++) {
                for(k=0; k<nw3; k++) fprintf(fp_write,"%lf\t%lf\t%lf\t%lf\t%lf\n",0.0,signal[0][l*nw3+k],signal[1][l*nw3+k]*fs_to_cm_m1/(2.*M_PI),signal[2][l*nw3+k],signal[3][l*nw3+k]);
            }
            
            fclose(fp_write);
            
            sprintf(nome_file_out,"FT_matrix_%s",nome_file_in);
            fp_write=fopen(nome_file_out,"w");
            
            for(k=0; k<nw3; k++) {
                for (l=0; l<nt2; l++) fprintf(fp_write,"%lf\t",signal[3][l*nw3+k]);
                fprintf(fp_write,"\n");
            }
            
            
            fclose(fp_write);

            
            printf("!! FILE WRITTEN !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            //------------------------------ PLOT RSPF AND FT ------------------------------//
            
            gnuplotPipe=popen(" gnuplot -persist", "w");
            fprintf(gnuplotPipe, "set xrange[%lf:%lf]\n",signal[0][0],signal[0][nw3*nt2-1]);
            fprintf(gnuplotPipe, "set yrange[%lf:%lf]\n",signal[1][0]*fs_to_cm_m1/(2.*M_PI),signal[1][nw3-1]*fs_to_cm_m1/(2.*M_PI));
            fprintf(gnuplotPipe, "set xlabel 't_2 (fs)'\n");
            fprintf(gnuplotPipe, "set ylabel 'w_3 (cm^{-1})'\n");
            fprintf(gnuplotPipe, "plot '%s' matrix u ($1*%lf):($2*%lf+%lf):3 w image title '' \n",nome_file_out,Dt2,Dw3*fs_to_cm_m1/(2.*M_PI),w3_min*fs_to_cm_m1/(2.*M_PI));
            fflush(gnuplotPipe);
            pclose(gnuplotPipe);
            
            
            gnuplotPipe=popen(" gnuplot -persist", "w");
            fprintf(gnuplotPipe, "plot '-' w l title 'real', '-' w l title 'imag' \n");
            l=0;
            for (i=0; i<nt3-1; i++) {
                fprintf(gnuplotPipe, "%lf\n",data_in[3][i+l*nt3]);
            }
            fprintf(gnuplotPipe, "e\n");
            for (i=0; i<nt3-1; i++) {
                fprintf(gnuplotPipe, "%lf\n",data_in[4][i+l*nt3]);
            }
            fprintf(gnuplotPipe, "e\n");
            fflush(gnuplotPipe);
            pclose(gnuplotPipe);

            //-------------------------------------------------------------------------------/*/
            
            
            //FREE SECTION
            for (i=0; i<5; i++) free(data_in[i]);
            for (i=0; i<4; i++) free(signal[i]);
            free(data_in);
            free(signal);
            
        }

        
        //FT OF 2DPP DATA
        if(action==3){
            
            int nt1, nt2, nw1, nt3, nw3;
            int ipn;
            int lt2_min, lt2_max;
            int choice;
            double sign;
            double Dt1, Dt2, Dt3, Dw1, Dw3, w3_min, w3_max, w1_min, w1_max;
            double **signal_auxiliary[2];
            double **signal[2];
            
            
            //--------------------------------- SELECT KI/KII/KIII/PP ---------------------------------//
            //you can easily compute also KIII     // +K1 +K2 -K3 -> how to visualize KIII? FT along t2 also?
            printf("KI, KII, KIII, PP (1, 2, 3, 0)?: ");
            scanf("%d",&choice);
            
            if (choice==1) printf("\nComputing KI responce (ks = -k1 +k2 +k3 & mirror) ...\n");
            else if (choice==2) printf("\nComputing KII responce (ks = +k1 -k2 +k3 & mirror)...\n");
            else if (choice==3) printf("\nComputing KIII responce (ks = +k1 +k2 -k3 & mirror)...\n");
            else printf("\nComputing PP (KI+KII) responce...\n");
            
            //modify!!!! -> print directly from pldm dynamics the file PP //
            //if(!choice) system("paste pldm.resp-k1 pldm.resp-k2 | awk \'{printf \"%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n\", $1, $2, $3, $4+$8,$5+$9}\' > pldm.resp-PP");
            //-----------------------------------------------------------------------------------------//
            
            
            
            
            printf("Insert RF file name (expected format t1 (fs), t2 (fs), t3 (fs), RE_f, IM_f): \n");
            if (choice==1) sprintf(nome_file_in,"pldm.resp-K1-++");
            else if (choice==2) sprintf(nome_file_in,"pldm.resp-K2+-+");
            else if (choice==3) sprintf(nome_file_in,"pldm.resp-K3");
            else sprintf(nome_file_in,"pldm.resp-PP");
            //scanf("%s",nome_file_in);
            
            
            //---------------------------------READ THE FILE---------------------------------//
            length=read_data_file(nome_file_in);
            
            data_in=malloc(5*sizeof(double*));
            for (i=0; i<5; i++) data_in[i]=malloc(length*sizeof(double));
            
            fp_read=fopen(nome_file_in,"r");
            
            nt2=1;
            nt3=0;
            for (i=0; i<length; i++) {
                fscanf(fp_read,"%lf\t%lf\t%lf\t%lf\t%lf\n",&data_in[0][i],&data_in[1][i],&data_in[2][i],&data_in[3][i],&data_in[4][i]);
                
                if((data_in[0][i]==data_in[0][0])&&(data_in[1][i]==data_in[1][0])) nt3++;
                if((data_in[1][i]!=data_in[1][0])&&(data_in[0][i]==data_in[0][0])&&(data_in[2][i]==data_in[2][0])) nt2++;
            }
            
            fclose(fp_read);
            
            nt1=length/nt3/nt2;
            Dt1=data_in[0][nt3]-data_in[0][0];
            Dt2=data_in[1][nt1*nt3]-data_in[1][0];
            Dt3=data_in[2][1]-data_in[2][0];
            
            printf("FILE METRICS:\n");
            printf("nt1: %.3f\tnt2: %.3f\tnt3: %.3f\n", (float)nt1,(float)nt2, (float)nt3); // %f for formatting reasons
            printf("Dt1: %.3f\tDt2: %.3f\tDt3: %.3f\n\n", Dt1 ,Dt2, Dt3);
            lt2_min=0;
            lt2_max=nt2;
            
            if(nt2>1) {
                printf("Specify the t2 interval: (0 -> %.2f; %d -> %.2f)\n",data_in[1][0],nt2-1,Dt2*nt2);
                scanf("%d\n%d\n",&lt2_min,&lt2_max); //maybe does not work properly
                if(lt2_min<0||lt2_max>=nt2) {
                    lt2_min=0;
                    lt2_max=nt2-1;
                }
                lt2_max+=1;
            }
            
            printf("!! READ FILES !!\n");

            //-------------------------------------------------------------------------------//

            
            //---------------------------------PARAMETERS---------------------------------//
            
            //AUTOMATIC w domain
            //freq in 1/fs
            w1_max=2.*M_PI/Dt1/2.; // x/2 because of nyquist
            w1_min=2.*M_PI/(Dt1*nt1)*2.; // *2
            nw1=nt1;
            
            w3_max=2.*M_PI/Dt3/2.; // x/2 because of nyquist
            w3_min=2.*M_PI/(Dt3*nt3)*2.; // *2
            nw3=nt3;
            
            //SET w domain
            printf("SET w1 DOMAIN cm^-1 (min, max, #): \n");
            scanf("%lf\t%lf\t%d", &w1_min, &w1_max, &nw1);
            w1_min=2.*M_PI*w1_min/fs_to_cm_m1;
            w1_max=2.*M_PI*w1_max/fs_to_cm_m1;
            
            printf("SET w3 DOMAIN cm^-1 (min, max, #): \n");
            scanf("%lf\t%lf\t%d", &w3_min, &w3_max, &nw3);
            w3_min=2.*M_PI*w3_min/fs_to_cm_m1;
            w3_max=2.*M_PI*w3_max/fs_to_cm_m1;
            
            Dw1=(w1_max-w1_min)/(nw1-1);
            Dw3=(w3_max-w3_min)/(nw3-1);
            
            //-------------------------------------------------------------------------------//

            
            //--------------------------------- COMPUTE FT ---------------------------------//
            //TRAPEZOID RULE! for every t2 and every w3 integrate over t3 (0-to-EOF), (RE_f*cos-IM_f*sin,IM_f*cos+RE_f*sin)
            
            
            for (ipn=0; ipn<2; ipn++) {
            
                //ALLOCATE MEMORY FOR THE signals//
                signal_auxiliary[ipn]=malloc(5*sizeof(double*));
                signal[ipn]=malloc(5*sizeof(double*));
                
                for (i=0; i<5; i++) {
                    signal[ipn][i]=calloc(nw1*nw3*nt2,sizeof(double));
                    signal_auxiliary[ipn][i]=calloc(nt1*nw3*nt2,sizeof(double));
                }
                ///////////////////////////////////
        
                //DURING SECOND CICLE READ THE OTHER FILE//
                if (ipn>0) {
                    printf("\nReading second file\n");
                    if (choice==1) sprintf(nome_file_in,"pldm.resp-K1+--");
                    else if (choice==2) sprintf(nome_file_in,"pldm.resp-K2-+-");
                    else if (choice==3) sprintf(nome_file_in,"pldm.resp-K3");
                    else sprintf(nome_file_in,"pldm.resp-PP");

                    fp_read=fopen(nome_file_in,"r");
                    
                    for (i=0; i<length; i++) fscanf(fp_read,"%lf\t%lf\t%lf\t%lf\t%lf\n",&data_in[0][i],&data_in[1][i],&data_in[2][i],&data_in[3][i],&data_in[4][i]);
                    
                    fclose(fp_read);
                }
                
                
                //////////////////// POSSIBLE MODIFICATION OF THE DATA ///////////////////
                for (i=0; i<length; i++) {
                    
                    if(i==0) printf("INTRODUCE 'COS' CORRECTION\n\n");
                    data_in[3][i]*=cos(M_PI*data_in[2][i]/(2.*data_in[2][nt3-1]))*cos(M_PI*data_in[0][i]/(2.*data_in[0][nt3*nt1-1]));
                    data_in[4][i]*=cos(M_PI*data_in[2][i]/(2.*data_in[2][nt3-1]))*cos(M_PI*data_in[0][i]/(2.*data_in[0][nt3*nt1-1]));
                    
                    
                }
                //////////////////////////////////////////////////////////////////////////
                
                
                
                //cicle over selected t2 time intervals
                for (l=lt2_min; l<lt2_max; l++) {
                    
                    printf("\nt2 = %.2f\n",data_in[1][l*nt1*nt3]);
                    
                    sign=(1.-2.*ipn);
                
                    //cicle over t1 times
                    for (j=0; j<nt1; j++) {
                        
                        //cicle over w3 frequencies
                        for(k=0; k<nw3; k++){
                            
                            //t1 time
                            signal_auxiliary[ipn][0][l*nw3*nt1+j*nw3+k]=data_in[0][l*nt3*nt1+j*nt3];
                            //t2 time
                            signal_auxiliary[ipn][1][l*nw3*nt1+j*nw3+k]=data_in[1][l*nt1*nt3];
                            //w3 freq
                            signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]=(w3_min+Dw3*k)*sign;
                            
                            //REAL
                            signal_auxiliary[ipn][3][l*nw3*nt1+j*nw3+k]=data_in[3][l*nt3*nt1+j*nt3]-data_in[4][l*nt3*nt1+(j+1)*nt3-1]*sin(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*(nt3-1))+data_in[3][l*nt3*nt1+(j+1)*nt3-1]*cos(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*(nt3-1));
                            signal_auxiliary[ipn][3][l*nw3*nt1+j*nw3+k]/=2;
                            //IMAGINARY
                            signal_auxiliary[ipn][4][l*nw3*nt1+j*nw3+k]=data_in[4][l*nt3*nt1+j*nt3]+data_in[3][l*nt3*nt1+(j+1)*nt3-1]*sin(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*(nt3-1))+data_in[4][l*nt3*nt1+(j+1)*nt3-1]*cos(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*(nt3-1));
                            signal_auxiliary[ipn][4][l*nw3*nt1+j*nw3+k]/=2;
                            
                            //cicle over the time t3 (exept 0 and last point, already accounted above)
                            for(i=1; i<nt3-1; i++){
                                
                                //REAL
                                signal_auxiliary[ipn][3][l*nw3*nt1+j*nw3+k]+=-data_in[4][i+l*nt3*nt1+j*nt3]*sin(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*i)+data_in[3][i+l*nt3*nt1+j*nt3]*cos(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*i);
                                //IMAGINARY
                                signal_auxiliary[ipn][4][l*nw3*nt1+j*nw3+k]+=data_in[3][i+l*nt3*nt1+j*nt3]*sin(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*i)+data_in[4][i+l*nt3*nt1+j*nt3]*cos(signal_auxiliary[ipn][2][l*nw3*nt1+j*nw3+k]*Dt3*i);
                                
                            }
                            
                            printf("\rCompleted %.2f\%%",((double)(j*nw3+k)+1.)*100./nt1/nw3);
                        }
                        
                    }
                    
                    printf("\n!! FT W3 COMPLETED !!\n");
                    
                    //--------------------------------------------------------//
                    
                    //KII
                    if (choice==2) sign=(1.-2.*ipn);
                    //KI/KIII -> modifica, se KIII la FT la devo fare su t2
                    else if (choice==1||choice==3) sign=-(1.-2.*ipn);

                    //cicle over w3 frequencies
                    for (j=0; j<nw3; j++) {
                        
                        //cicle over w1 frequencies
                        for(k=0; k<nw1; k++){
                            
                            //w1 freq
                            signal[ipn][0][l*nw3*nw1+k*nw3+j]=(w1_min+Dw1*k)*sign;
                            //t2 time
                            signal[ipn][1][l*nw3*nw1+k*nw3+j]=signal_auxiliary[ipn][1][l*nt1*nw3];
                            //w3 freq
                            signal[ipn][2][l*nw3*nw1+k*nw3+j]=signal_auxiliary[ipn][2][l*nt1*nw3+j];
                            
                            
                            //REAL
                            signal[ipn][3][l*nw3*nw1+k*nw3+j]=signal_auxiliary[ipn][3][l*nw3*nt1+j]-signal_auxiliary[ipn][4][l*nw3*nt1+(nt1-1)*nw3+j]*sin(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*(nt1-1))+signal_auxiliary[ipn][3][l*nw3*nt1+(nt1-1)*nw3+j]*cos(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*(nt1-1));
                            signal[ipn][3][l*nw3*nw1+k*nw3+j]/=2;
                            //IMAGINARY
                            signal[ipn][4][l*nw3*nw1+k*nw3+j]=signal_auxiliary[ipn][4][l*nw3*nt1+j]+signal_auxiliary[ipn][3][l*nw3*nt1+(nt1-1)*nw3+j]*sin(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*(nt1-1))+signal_auxiliary[ipn][4][l*nw3*nt1+(nt1-1)*nw3+j]*cos(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*(nt1-1));
                            signal[ipn][4][l*nw3*nw1+k*nw3+j]/=2;
                            
                            //cicle over the time t1 (exept 0 and last point, already accounted above)
                            for(i=1; i<nt1-1; i++){
                                
                                //REAL
                                signal[ipn][3][l*nw3*nw1+k*nw3+j]+=-signal_auxiliary[ipn][4][i*nw3+l*nw3*nt1+j]*sin(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*i)+signal_auxiliary[ipn][3][i*nw3+l*nw3*nt1+j]*cos(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*i);
                                //IMAGINARY
                                signal[ipn][4][l*nw3*nw1+k*nw3+j]+=signal_auxiliary[ipn][3][i*nw3+l*nw3*nt1+j]*sin(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*i)+signal_auxiliary[ipn][4][i*nw3+l*nw3*nt1+j]*cos(signal[ipn][0][l*nw3*nw1+k*nw3+j]*Dt1*i);
                                
                            }
                            
                            printf("\rCompleted %.2f\%%",((double)(j*nw1+k)+1.)*100./nw1/nw3);
                            
                        }

                    }
                    
                    printf("\n!! FT W1 COMPLETED !!\n");
                    
                }
                
                printf("\n!! COMPLETE FT PERFORMED !!\n\n");
                
                
            }
            
            //change name to have same output
            if (choice==1) sprintf(nome_file_in,"pldm.resp-K1");
            else if (choice==2) sprintf(nome_file_in,"pldm.resp-K2");
            else if (choice==3) sprintf(nome_file_in,"pldm.resp-K3");
            else sprintf(nome_file_in,"pldm.resp-PP");

            
            //-------------------------------------------------------------------------------//
            
            
            //--------------------------------- SAVE FILE ---------------------------------//
            
            //IF NEEDED SHIFT THE FREQUENCY DOMAIN (E.G. IN SPECTRON THE CARRYING FREQUENCY IS REMOVED IN TIME DOMAIN AND REINTRODUCED AS A SHIFTING IN THE FREQ DOMAIN)
            
            
            ////FREE AUXILIARY FUNCTION////
            for (ipn=0; ipn<2; ipn++) {
                for (i=0; i<5; i++) free(signal_auxiliary[ipn][i]);
                free(signal_auxiliary[ipn]);
            }
            /////////////////////////
                
            
            //KI note -signal[0][0][j*nw3+k+l*nw3*nw1]
            if (choice==1) {
           
                sprintf(nome_file_out,"FT_%s",nome_file_in);
                fp_write=fopen(nome_file_out,"w");
                
                for (l=lt2_min; l<lt2_max; l++) {
                    for (j=0; j<nw1; j++){
                        for(k=0; k<nw3; k++) {
                            fprintf(fp_write,"%lf\t%lf\t%lf\t%lf\t%lf\n",-signal[0][0][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[0][1][j*nw3+k+l*nw3*nw1],signal[0][2][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[0][3][j*nw3+k+l*nw3*nw1]+signal[1][3][j*nw3+k+l*nw3*nw1],signal[0][4][j*nw3+k+l*nw3*nw1]-signal[1][4][j*nw3+k+l*nw3*nw1]);
                        }
                    }
                }
                
                fclose(fp_write);
                
            }
            
            //KII note +signal_pos[0][0][j*nw3+k+l*nw3*nw1]
            if (choice==2) {
                
                sprintf(nome_file_out,"FT_%s",nome_file_in);
                fp_write=fopen(nome_file_out,"w");
                
                for (l=lt2_min; l<lt2_max; l++) {
                    for (j=0; j<nw1; j++){
                        for(k=0; k<nw3; k++) {
                            fprintf(fp_write,"%lf\t%lf\t%lf\t%lf\t%lf\n",signal[0][0][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[0][1][j*nw3+k+l*nw3*nw1],signal[0][2][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[0][3][j*nw3+k+l*nw3*nw1]+signal[1][3][j*nw3+k+l*nw3*nw1],signal[0][4][j*nw3+k+l*nw3*nw1]-signal[1][4][j*nw3+k+l*nw3*nw1]);
                        }
                    }
                }
                
                fclose(fp_write);
                
                
            }

            //KIII note -signal[0][2][j*nw3+k+l*nw3*nw1] e inverted signal
            if (choice==3) {
                
                sprintf(nome_file_out,"FT_%s",nome_file_in);
                fp_write=fopen(nome_file_out,"w");
                
                for (l=lt2_min; l<lt2_max; l++) {
                    for (j=0; j<nw1; j++){
                        for(k=0; k<nw3; k++) {
                            fprintf(fp_write,"%lf\t%lf\t%lf\t%lf\t%lf\n",-signal[0][0][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[0][1][j*nw3+k+l*nw3*nw1],signal[0][2][j*nw3+k+l*nw3*nw1]*fs_to_cm_m1/(2.*M_PI),signal[1][3][j*nw3+k+l*nw3*nw1]+signal[0][3][j*nw3+k+l*nw3*nw1],signal[1][4][j*nw3+k+l*nw3*nw1]-signal[0][4][j*nw3+k+l*nw3*nw1]);
                        }
                    }
                }
                
                fclose(fp_write);
                
                
            }

            
            
            printf("!! FILE WRITTEN !!\n");
            
            //-------------------------------------------------------------------------------//
            
            
            //------------------------------ PLOT RSPF AND FT ------------------------------//
            //CHANGE IN CASE OF SEVERAL T2 TIMES
            /*
            gnuplotPipe=popen(" gnuplot -persist", "w");
            fprintf(gnuplotPipe, "set xlabel 'w_1 (cm^{-1})'\n");
            fprintf(gnuplotPipe, "set ylabel 'w_3 (cm^{-1})'\n");
            for (l=0; l<nt2; l++) {
                fprintf(gnuplotPipe, "set title 't2 = %.3f'\n",signal[1][l*nw3*nw1]);
                fprintf(gnuplotPipe, "plot 'FT_PP_pldm.resp' every 1:1:%d:0:%d:0 u 1:3:(-$5) w image title '' \n",l*nw1*nw3,(l+1)*nw1*nw3-1);
                fflush(gnuplotPipe);
                poll(0, 0, 3000);
            }
            
            pclose(gnuplotPipe);
             */
            
            //-------------------------------------------------------------------------------/*/
            
            
            //FREE SECTION
            for (i=0; i<5; i++) free(data_in[i]);
            for (i=0; i<5; i++) {
                free(signal[0][i]);
                free(signal[1][i]);
            }
            free(data_in);
            
            free(signal[0]);
            free(signal[1]);
            
        }

        
        
    }//end of while
    
    
    return 0;
    
}



int read_data_file(char *nome_file){
    
    int i, length=0;
    int ch;
    FILE *fp_read;
    
    fp_read=fopen(nome_file,"r");
    
    while (!feof(fp_read)) {
        ch=fgetc(fp_read);
        if(ch == '\n') length++;
    }
    
    fclose(fp_read);
    
    return length;
}


/*
 
 for(k=0; k<n_w; k++){
 
 //REAL
 signal[1][k]=data_in[1][0]-data_in[2][length-1]*sin(2.*(1./length)*M_PI*k*(length-1))+data_in[1][length-1]*cos(2.*(1./length)*M_PI*k*(length-1));
 signal[1][k]/=2;
 //IMAGINARY
 signal[2][k]=data_in[2][0]+data_in[1][length-1]*sin(2.*(1./length)*M_PI*k*(length-1))+data_in[2][length-1]*cos(2.*(1./length)*M_PI*k*(length-1));
 signal[2][k]/=2;
 
 //cicle over the time
 for(i=1; i<length-1; i++){
 
 //REAL
 signal[1][k]+=-data_in[2][i]*sin(2.*(1./length)*M_PI*k*i)+data_in[1][i]*cos(2.*(1./length)*M_PI*k*i);
 //IMAGINARY
 signal[2][k]+=data_in[1][i]*sin(2.*(1./length)*M_PI*k*i)+data_in[2][i]*cos(2.*(1./length)*M_PI*k*i);
 
 }
 
 printf("\rCompleted %.2f\%%",((double)k+1)*100./length);
 
 }

 
 */
