#include <stdio.h>
#include <math.h>

int main(){
double norm_T=800*0.472/0.56;
int norm_x=2;
double norm_t=pow(0.02,2)*3686*1081/0.56;
int DIM=100;
double T_amb=2*0.56*(36.5+273.13)/(pow(40,2)*0.472);;
double delta_x=1./(DIM-1);
//double delta_t=delta_x;
//double delta_t=delta_x*0.5;
double delta_t=pow(delta_x,2)*0.25;
int T_it=1000;
double T[T_it][DIM];
double E[T_it][DIM];
double gamma=delta_t/pow(delta_x,2);
double x=0;
int temps_aturat=0;
char buffer_0[150];
FILE *fitxer;
sprintf(buffer_0,"mz%d.txt",0);
fitxer=fopen(buffer_0,"w+");
for(int i=0;i<DIM+1;i++){//S'imposen les condicions inicials i es crea el primer fitxer
    T[0][i]=T_amb;
    fprintf(fitxer, "%.30lf\t %lf \n",T[0][i], i* delta_x);
}
fclose(fitxer);
for(int k=1;k<=T_it;k++){//Es crea el bucle temporal
    char buffer[150];//Es crea la resta de fitxers fins la l�nia 34
    FILE *fitxer;
    sprintf(buffer,"mz%d.txt",k);
    fitxer=fopen(buffer,"w+");
    fprintf(fitxer, "%lf \n", norm_t*k*delta_t);
    double T_anterior[DIM];
    for(int l=0;l<DIM;l++){//Imposen que la temperatura amb la qual iterarem ser� aproximadament igual que la temperatura
    //anterior per comen�ar la nostre aproximaci� de Guass-Seidel
        T[k][l]=T[k-1][l];
        T_anterior[l]=T[k][l];
    }
double h=10;//Es crea una varaible que ens aturar� gauss seidel quan sigui prou petita
    double Dif[DIM];//Aquest vector ser� la difer�ncia entre el pas actual i el pas anterior de Gauss-Seidel
    while(h>0.0000000001){// Es crea el bucle de Gauss-Seidel, que imposarem que pari quan la difer�ncia
            //m�xima entre cada iteraci� sigui prou petita
        for(int j=1;j<DIM-1;j++){//S'estableix la equaci� que s'ha trobat per solucionar el problema
            T[k][j] = T[k-1][j] + (gamma/2)*(T[k-1][j+1]-2*T[k-1][j]+T[k-1][j-1] + T[k][j+1]-2*T[k][j]+T[k][j-1]) + delta_t;
        }
        for(int j=0;j<DIM;j++){//Es mira la difer�ncia entre la  iteraci� actual i la anterior
            Dif[j]=fabs(T[k][j]-T_anterior[j]);
            T_anterior[j]=T[k][j];
        }
        double max=Dif[0];
        for(int j=1;j<DIM;j++){//Es mira el m�xim de les difer�ncies entre els dos passos de Gauss-Seidel
            if(Dif[j]>Dif[j-1]){
                max=Dif[j];
            }
        h=max;//Si aquest m�xim �s menor al n�mero que hi ha en el bucle, aquest s'aturar�
        }
    }
    for(int a=0;a<DIM;a++){
        fprintf(fitxer, "%.30lf\t %lf \n",norm_T*T[k][a]-273.13, norm_x*a* delta_x);
    }
    fclose(fitxer);
    for(int j=1;j<DIM;j++){//creo un nou bucle per imposar les condicions l�mit del problema
        if(j*delta_x<0.375 || j*delta_x>0.625){//Agafo nom�s els punts de la zona que no est� sana
            double a=norm_T*T[k][j]-273.13;
            if(a>=(50)){//creo un if que em surt del bucle temporal si la temperatura de la regi� sana pasa de 50 graus
                goto res;//la funci� "goto" envia fora del bucle temporal
            }
        }
    }
    x=0;
    temps_aturat+=1;
}
res: printf("La temperatura es major a 80 graus i el temps en que s'ha parat %cs %lf i el fitxer %cltim sera %i",130,norm_t*(temps_aturat)*delta_t,163,temps_aturat);
//amb aquest �ltim pas ja hem sortit del bucle i ens ha dit quin ser� el nostre �ltim fitxer
return 0;
}
