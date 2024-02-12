#include <stdio.h>
#include <math.h>



int main(){
int DIM=100;
double norm_T=800*0.472/0.56;
int norm_x=2;
double norm_t=pow(0.02,2)*3686*1081/0.56;
double T_amb=2*0.56*(36.5+273.13)/(pow(40,2)*0.472);
double delta_x=1./(DIM-1);
double delta_t=pow(delta_x,2)*0.25;
//double delta_t=pow(delta_x,2)*0.49;// són les diferents disctretitzacions del temps
//double delta_t=pow(delta_x,2)*0.51;
int t=1000;
double T[t][DIM];
int temps_aturat=0;
double x=0;
double gamma=delta_t/pow(delta_x,2);
char buffer_0[150];
FILE *fitxer;
sprintf(buffer_0,"mz%d.txt",0);
fitxer=fopen(buffer_0,"w+");
for(int i=0;i<DIM+1;i++){//S'imposen les condicions inicials i es crea el primer fitxer
    T[0][i]=T_amb;
    fprintf(fitxer, "%.30lf\t %lf \n",norm_T*T[0][i]-273.13, norm_x*i* delta_x);
}
fclose(fitxer);
for(int k=0;k<=t;k++){//aquest serà el bucle que iterarà sobre el temps
    char buffer[150];//des d'aquí fins la línia 35 només creo els fitxers en els quals em guardaré les dades
    FILE *fitxer;
    sprintf(buffer,"mz%d.txt",k+1);
    fitxer=fopen(buffer,"w+");
    fprintf(fitxer, "%lf \n", norm_t*(k+1)*delta_t);
    for(int j=0;j<=DIM-1;j++){//aquí ja començo amb el bucle espaial, que a dins tindrà un if per establir les condicions de contorn
        if(j==0){//estableixo que pel primer punt de cada temps tingui la temperatura ambient
            T[k+1][j]=T_amb;
            fprintf(fitxer, "%5.4lf \t %lf \n", norm_T*T[k+1][j]-273.13, norm_x*x);
            x+=delta_x;
        }
        else if(j==DIM-1){//estableixo que per l'últim punt de cada temps tingui la temperatura ambient
            T[k+1][j]=T_amb;
            fprintf(fitxer, "%5.4lf \t %lf \n", norm_T*T[k+1][j]-273.13, norm_x*x);
            x+=delta_x;
        }
        else{//pels altres punts poso la equació que he trobat amb el mètode d'euler explícit
            T[k+1][j]=gamma*(T[k][j+1]+T[k][j-1])+T[k][j]*(1-2*gamma)+delta_t;
            fprintf(fitxer, "%5.4lf \t %lf \n", norm_T*T[k+1][j]-273.13, norm_x*x);
            x+=delta_x;
        }
    }
    x=0;
    fclose(fitxer);
    for(int j=0;j<DIM;j++){//creo un nou bucle per imposar les condicions límit del problema
        if(j*delta_x<0.75/2 || j*delta_x>1.25/2){//agafo només els punts de la zona que no està sana
            double a=norm_T*T[k][j]-273.13;
            //printf("%lf\n",a);
            if(a>50){//creo un if que em surt del bucle temporal si la temperatura de la regió sana pasa de 50 graus
                printf("LA temperatura es major a 50 graus en el temps %lf\n", delta_t*k*norm_t);
                goto res;//la funció "goto" m'envia fora del bucle temporal
            }
        }


    }
temps_aturat+=1;//aquesta variable després ens tornarà en quin temps s'ha resolt el problema
}
res: printf("LA temperatura es major a 50 graus i el temps en que s'ha parat %cs %lf i el fitxer %cltim sera %i",130,norm_t*(temps_aturat-1)*delta_t,163,temps_aturat-1);
//amb aquest últim pas ja hem sortit del bucle i ens ha dit quin serà el nostre últim fitxer
return 0;
}
