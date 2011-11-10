#include "Matrix.h"
#include <iostream>
#include <cmath>

// main()-Funktion, no parameter
int main(int argc, char* argv[]){


   //anzahl trainingselemente
  int P = 11;
  //komplaexitaet
  int M = 5;
  
  double x[P];
  double t[P];

  //trainingset berechnen
  for(int i=0;i<P;i++){
          x[i]=(((2*i-1)*M_PI)/P);
          t[i]=cos(x[i]/2);
          }

  
  

  Matrix MC(M+1, M+2);
  MC.initialize(0.0);

 //a_{km}
   for(int row=0;row<=M;row++){
          for(int column=0;column<=M;column++){
                float temp=0;  
                          for(int j=0;j<P;j++)
                 temp+= pow(x[j],row+column);
                  
          double* matrixRow=MC[row];
          matrixRow[column]=temp;
                  }
          }          
  
  
   
  
  
  //b_k
  for(int i=0;i<=M;i++){
          float temp=0;
          for(int j=0;j<P;j++)
                 temp+= t[j]*pow(x[j],i);
          
          double* row=MC[i];
          row[M+1]=temp;
    }
    
std::cout << MC;
//loesen
   for(int row=0;row<=M;row++){

          double* matrixRow=MC[row];
          double divisor= matrixRow[row];

          //ueber die reihe iterieren und alles durch die zahl der aktuellen pos teilen
          for(int position=0;position<=M+1;position++)
                  matrixRow[position]=matrixRow[position]/divisor;

          //diese Reihe von allen anderen subtrahieren.
          for(int rows=0;rows<=M;rows++){
                  //die zeile die von allen anderen subtrahiert wird ueberspringen
                  if(rows==row)
                        continue;
                  
                  double* currentRow=MC[rows];     
                  double unterschied=currentRow[row]; 
                  
                  //jedes element aus der einen zeile von der anderen abziehen unter 
                  //beruecksichtigung des unterschiedes zwischen dem 1er element und dem 
                  //korrespondierenden in der anderen zeile 
                  //(damit die position in der in einer zeile 1 steht in der 
                  //anderen 0 wird)
                  for(int column=0;column<=M+1;column++){
                       currentRow[column]-=unterschied*matrixRow[column];   
                  }
          }        

    }    

std::cout << MC;

  return 0;
}
