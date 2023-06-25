#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Inicialização para o modelo de água SPC

void main(){
    
    //Variáveis
    
    // MOLECULAS
    int N = 150; // total de particulas
    int bonds = 100;
    int angles = 50;
    int a = 1; // atom 1 = hy ; atom 2 = ox
    float q[2] = {0.4238, -0.8476};
    float mol[3][4];

    // CAIXA DE SIMULACAO
    double x = 0.1, y= 0.1, z = 0.1; // coordenadas das particulas
    double lx = 40, ly = 40, lz = 40; // dimensoes da caixa
    double lx0 = -5, ly0 = -5, lz0 = -5; // dimensoes da caixa

    //double lx = 44.165558, ly = 47.066074, lz = 44.193867; // dimensoes da caixa
    //double lx0 = -0.75694412, ly0 = 0.38127473, lz0 = 0.17900842; // dimensoes da caixa


    double rho = N/(lx*ly*lz);
    double rho_l = cbrt(rho);
    int Nx = round(rho_l*lx);
    int Ny = round(rho_l*ly);
    int Nz = round(rho_l*lz);
    int n3 = ceil(cbrt(N))-1;
    double espaco = 1/rho_l;

    float dist = 1;
    int i = 0, j = 0, k = 0;

    int mol_counter = 1;
    int r = 1;
    
    
    // Arquivos para impressão
    FILE *arq_Init;
    arq_Init = fopen("init.syst", "w");
    
    fprintf(arq_Init, "#LAMMPS SPC-WATER\n");
    fprintf(arq_Init, "        \n%d atoms        \n%d bonds        \n%d angles", N, bonds, angles);
    fprintf(arq_Init, "           \n2 atom types           \n1 bond types           \n1 angle types\n");
    
    fprintf(arq_Init, "0.0   %.1f   xlo xhi\n0.0   %.1f   ylo yhi\n0.0   %.1f   zlo zhi\n\n", lx, ly, lz);
    
    fprintf(arq_Init, "Masses\n\n     1  1.0079401\n     2  15.999400\n\n");
    
    fprintf(arq_Init, "Atoms\n\n");
    
    printf("Densidade Total: %f\n", rho);
    printf("Densidade Lateral: %f\n", rho_l);
    printf("Numero de particulas: x = %d y = %d z = %d sum = %d\n", Nx, Ny, Nz, Nx*Ny*Nz);
    printf("Distancias entre particulas: %f\n", espaco);


    // Cria 1 molécula
    a = 2;
    for(int n = 0 ;  n <= 2 ; n++){

        mol[n][0] = a; //guarda o tipo de partícula

        y = n*dist + 0.1;
        
        if(n == 2){
            y = dist*cos(1.910632); // angulo entre moleculas é 109.4712
            z = dist*sin(1.910632);
        }

        // guarda as posições
        mol[n][1] = x;
        mol[n][2] = y;
        mol[n][3] = z;

        a = 1;
    }

    int p = 0;
    for(int n = 1 ; n <= N ; n++){

        x = mol[p][1] + i*(dist + espaco);
        y = mol[p][2] + j*(dist + espaco);
        z = mol[p][3] + k*(dist + espaco);

        fprintf(arq_Init, "         %3d           %d           %d   %10.4f        %10.16f        %10.16f        %10.16f     \n", n, r, (int) mol[p][0], q[((int) mol[p][0]) - 1], x, y, z);

        p++;
        if(p == 3){
            p = 0;
            r++;

            i++;
            if(i == n3){
                i = 0;
                j++;
                if(j == n3){
                    j = 0;
                    k++;
                }
            }
        }
    }


    fprintf(arq_Init, "\nBonds\n\n");
    int b = 1;


    for(int a = 1 ; a <= bonds; a = a + 2){

        fprintf(arq_Init, "         %4d        1        %4d        %4d\n", a, b, b+1);
        fprintf(arq_Init, "         %4d        1        %4d        %4d\n", a+1, b, b+2);
        b = b + 3;
    }

    fprintf(arq_Init, "\nAngles\n\n");
    b = 1;

    for(int a = 1 ; a <= angles; a++){

        fprintf(arq_Init, "         %4d        1        %4d        %4d        %4d\n", a, b+1, b, b+2);
        b = b + 3;
    }

    fclose(arq_Init);
}

    /*
    for(int n = 1 ; n <= N ; n++){
        
        x = i*espaco + 0.1;
        y = j*espaco + 0.1;
        z = k*espaco + 0.1;
        
        fprintf(arq_Init, "         %3d           %d           %d   %10.4f        %10.16f        %10.16f        %10.16f     \n", n, r, a, q, x, y, z);

        i++;
        if(i == n3){
            i = 0;
            j++;
            if(j == n3){
                j = 0;
                k++;
            }
        }

        if(mol_counter == 2){
            // pois cada ligação tem 2 hidrogênios
            a = 2;
            q = -0.8476;
            mol_counter++;
        }
        else if(mol_counter == 3){
            a = 1;
            q = 0.4238;
            mol_counter = 1;
            r++;
        }
        else{
            mol_counter++;
        }

    // Cria 1 molécula
    float dist = 1;
    x = 0.1;
    z = 0.1;
    a = 2;
    for(int n = 0 ;  n <= 2 ; n++){

        y = n*dist + 0.1;

        if(n == 2){
            y = dist*cos(1.910632); // angulo entre moleculas é 109.4712
            z = dist*sin(1.910632);
        }

        fprintf(arq_Init, "         %3d           %d           %d   %10.4f        %10.16f        %10.16f        %10.16f     \n", n, r, a, q, x, y, z);
        a = 1;
    }
    */