#include <stdio.h>
#include <math.h>
int main(){
double norm_T=800*0.472/0.56;
double norm_t=pow(0.02,2)*3686*1081/0.56;
int norm_x=2;
double x_max=1,x_min=0;
double PI = 3.14159265358979323846;
int DIM=100,q=0;
double T_amb=2*0.56*(36.5+273.13)/(pow(40,2)*0.472);
double delta_x=1./(DIM-1);
//double delta_t=delta_x*0.5;
//double delta_t=delta_x;
double delta_t=pow(delta_x,2)*0.25;
double t_max=0.025, t_min=0,t=delta_t;
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
    E[0][i]=T[0][i]-T_amb;
    fprintf(fitxer, "%.30lf\t %lf \n",norm_T*E[0][i], norm_x*i* delta_x);
}
fclose(fitxer);

for(int k=1;k<=T_it;k++){//Es comença el bucle temporal
    char buffer[150];// Es creen els altres fitxers fins la línia 39
    FILE *fitxer;
    sprintf(buffer,"mz%d.txt",k);
    fitxer=fopen(buffer,"w+");
    fprintf(fitxer, "%lf \n", norm_t*k*delta_t);
    double T_anterior[DIM];
    for(int l=0;l<DIM+1;l++){//Imposen que la temperatura amb la qual iterarem serà aproximadament igual que la temperatura
    //anterior per començar la nostre aproximació de Guass-Seidel
        T[k][l]=T[k-1][l];
        T_anterior[l]=T[k][l];
    }
    double h=10;//Es crea una varaible que ens aturarà gauss seidel quan sigui prou petita
    double Dif[DIM];//Aquest vector serà la diferència entre el pas actual i el pas anterior de Gauss-Seidel
    while(h>0.0000000001){// Es crea el bucle de Gauss-Seidel, que imposarem que pari quan la diferència
            //màxima entre cada iteració sigui prou petita
        for(int j=1;j<DIM-1;j++){//S'estableix la equació que s'ha trobat per solucionar el problema
            T[k][j]=(gamma*(T[k][j-1]+T[k][j+1])+delta_t+T[k-1][j])/(1+2*gamma);
        }
        for(int j=0;j<DIM;j++){//Es mira la diferència entre la  iteració actual i la anterior
            Dif[j]=fabs(T[k][j]-T_anterior[j]);
            T_anterior[j]=T[k][j];
        }
        double max=Dif[0];
        for(int j=1;j<DIM;j++){//Es mira el màxim de les diferències entre els dos passos de Gauss-Seidel
            if(Dif[j]>Dif[j-1]){
                max=Dif[j];
            }
        h=max;//Si aquest màxim és menor al número que hi ha en el bucle, aquest s'aturarà
        }
    }
    for(int a=0;a<DIM;a++){//Es mira la diferència entra la solució analítica i la numèrica
        double sum=0;
        for(q=0;q<1000;q++){//Es crea la solució analítica per la posició i temps en les que ens trobem
            int f=2*q+1;
            sum+=(1-exp(-pow(f,2)*pow(PI,2)*t))*sin(f*PI*x)/(pow(PI,3)*pow(f,3));
        }
        double g=T_amb+4*sum;
        fprintf(fitxer, "%.30lf \t %lf \n", norm_T*fabs(T[k][a]-g),norm_x* a* delta_x);
        x+=delta_x;
    }
    x=0;
    t+=delta_t;
    fclose(fitxer);
     for(int j=1;j<DIM;j++){//creo un nou bucle per imposar les condicions límit del problema
        if(j*delta_x<0.375 || j*delta_x>0.625){//Agafo només els punts de la zona que no està sana
            double a=norm_T*T[k][j]-273.13;
            if(a>=(50)){//creo un if que em surt del bucle temporal si la temperatura de la regió sana pasa de 50 graus
                goto res;//la funció "goto" m'envia fora del bucle temporal
            }
        }
    }
    temps_aturat+=1;
}
res: printf("LA temperatura es major a 80 graus i el temps en que s'ha parat és %lf i el fitxer ultim sera %i",norm_t*(temps_aturat)*delta_t,temps_aturat);
//amb aquest últim pas ja hem sortit del bucle i ens ha dit quin serà el nostre últim fitxer
return 0;
}

