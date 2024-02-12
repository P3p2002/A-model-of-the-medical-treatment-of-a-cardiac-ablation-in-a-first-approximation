#include <stdio.h>
#include <math.h>

int main(){
double norm_T=800*0.472/0.56;
int norm_x=2;
double norm_t=pow(0.02,2)*3686*1081/0.56;
double x_max=1,x_min=0;
double PI = 3.14159265358979323846;
int DIM=100,q=0;
double T_amb=2*0.56*(36.5+273.13)/(pow(40,2)*0.472);
double delta_x=1./(DIM-1);
double delta_t=pow(delta_x,2)*0.25;
//double delta_t=pow(delta_x,2)*0.49;
double t_max=0.025, t_min=0,t=delta_t;
int T_it=1000;
double T[T_it][DIM];
double E[T_it][DIM];
double gamma=delta_t/pow(delta_x,2);
double x=0;
double temps_aturat;
char buffer_0[150];
FILE *fitxer;
sprintf(buffer_0,"mz%d.txt",0);
fitxer=fopen(buffer_0,"w+");
for(int i=0;i<DIM+1;i++){//S'imposen les condicions inicials i es crea el primer fitxer
    T[0][i]=T_amb;
    E[0][i]=fabs(T[0][i]-T_amb);
    fprintf(fitxer, "%.30lf\t %lf \n",norm_T*E[0][i], norm_x*i* delta_x);
}
fclose(fitxer);

for(int k=0;k<=T_it;k++){//començo el bucle temporal
    char buffer[150];//des d'aquí fins la línia 38 només creo els fitxers en els quals em guardaré les dades
    FILE *fitxer;
    sprintf(buffer,"mz%d.txt",k+1);
    fitxer=fopen(buffer,"w+");
    fprintf(fitxer, "%lf \n", norm_t*(k+1)*delta_t);
    for(int j=0;j<=DIM-1;j++){//començo el bucle espaial

        double sum=0;//des d'aqui fins la línia 47 trobo la solució analítica per el punt i temps dels pasos creats i espaials que hi ha en els fors
        for(q=0;q<1000;q++){
        int f=2*q+1;
        sum+=(1-exp(-pow(f,2)*pow(PI,2)*t))*sin(f*PI*x)/(pow(PI,3)*pow(f,3));
        }
        double f=T_amb+4*sum;

        if(j==0){//imposo la condició de contorn d'un extrem i miro l'error en aquest punt
            T[k+1][j]=T_amb;
            E[k+1][j]=T[k+1][j]-f;
            fprintf(fitxer, "%.10lf \t %lf \n", fabs(norm_T*E[k+1][j]),norm_x*x);
            x+=delta_x;
        }
        else if(j==DIM-1){//imposo la condició de contorn de l'altre extrem i miro l'error en aquest punt
            T[k+1][j]=T_amb;
            E[k+1][j]=T[k+1][j]-f;
            fprintf(fitxer, "%.10lf \t %lf \n", fabs(norm_T*E[k+1][j]),norm_x*x);
            x+=delta_x;
        }
        else{//pels altres punts miro quina solució tinc i també el seu errro
            T[k+1][j]=gamma*(T[k][j+1]+T[k][j-1])+T[k][j]*(1-2*gamma)+delta_t;
            E[k+1][j]=T[k+1][j]-f;
            fprintf(fitxer, "%.10lf \t %lf \n", fabs(norm_T*E[k+1][j]),norm_x*x);
            x+=delta_x;
        }
    }
    x=0;
    t+=delta_t;
    fclose(fitxer);
    for(int j=0;j<DIM;j++){//creo un nou bucle per imposar les condicions límit del problema
        if(j*delta_x<0.75/2 || j*delta_x>1.25/2){//agafo només els punts de la zona que no està sana
            double a=norm_T*T[k][j]-273.13;
            //printf("%lf\n",a);
            if(a>50){//creo un if que em surt del bucle temporal si la temperatura de la regió sana pasa de 50 graus
                goto res;//la funció goto m'envia fora del bucle temporal
            }
        }


    }
    temps_aturat+=1;
}
res: printf("LA temperatura es major a 50 graus i el temps en que s'ha parat és %lf i el fitxer ultim sera %lf", norm_t*(temps_aturat-1)*delta_t, temps_aturat-1);
//amb aquest últim pas ja hem sortit del bucle i ens ha dit quin serà el nostre últim fitxer
return 0;
}

