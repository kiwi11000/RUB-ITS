#include "Matrix.h"
#include <iostream>
#include <cmath>


double calcAndReturnError(double* x,double* t,double* x1,double* t1,int M,int P);

// main()-Funktion, no parameter
int main(int argc, char* argv[]){


   //anzahl trainingselemente
  int P = 11;
  //komplaexitaet
  int M = 5;
  int M_MAX=20;
  
  //trainingsset
  double x[P];
  double t[P];
  
  //testset
  double xDach[P];
  double tDach[P];


  //trainingset berechnen
  for(int i=0;i<P;i++){
          x[i]=(((2*i-1)*M_PI)/P);
          //x[i]=pow(M_E,pow(i+1,2)/(2*pow(0.6,2)));
          std::cout << x[i] << std::endl;
          t[i]=cos(x[i]/2);
          }

   //testset berechnen
     for(int i=0;i<P;i++){
          xDach[i]=((2*i*M_PI)/P)-M_PI;
          tDach[i]=cos(xDach[i]/2);
          }


   for(int i=0;i<M_MAX;i++){
     std::cout << std::endl << std::endl << "Komplaexitaet M= " << i << std::endl;

     //Fehler des Trainingssets berechnen (param x, t sind die Trainingssetdaten)  
     std::cout << "Error vom Trainingsset: " << calcAndReturnError(x,t, xDach, tDach, i, P)<< std::endl;
     //Fehler des Testsets berechnen (param xDach, tDach sind die Testsetdaten)
     //std::cout << "Error vom Testset: " <<calcAndReturnError(xDach,tDach, i, P)<<std::endl;
   }



  return 0;
}

//
double calcAndReturnError(double* x,double* t,double* x1,double* t1,int M,int P){
       
  Matrix MC(M+1, M+2);
  MC.initialize(0.0);

 //a_{km}
   for(int row=0;row<=M;row++){
          for(int column=0;column<=M;column++){
                float temp=0.0;  
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
          //std::cout << row[M+1] << " " << temp << std::endl;
    }
    
    
    //std::cout << MC;

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

    //std::cout << MC;

    //array fuer alle w werte(koennte man vielleicht auch via matrix klasse machen)
   double w[M+1];
   //in array w die w werte speichern
   for(int row=0;row<=M;row++){
        double* currentRow=MC[row];
        w[row]=currentRow[M+1];
        //std::cout << M << " " << " " << currentRow[M+1] << " " << w[row] << std::endl;
   }

   //E(w)
   double fehler=0.0;
       for(int p=0;p<P;p++){
       //y(w,x)  
       double y=0.0;
       for(int innerSum=0;innerSum<=M;innerSum++)
             y+=w[innerSum]*pow(x[p],innerSum);
     
       fehler+=pow(y-t[p],2);
       }

std::cout << fehler << std::endl;
   //E(w)
   fehler=0.0;
       for(int p=0;p<P;p++){
       //y(w,x)  
       double y=0.0;
       for(int innerSum=0;innerSum<=M;innerSum++)
             y+=w[innerSum]*pow(x1[p],innerSum);
     
       fehler+=pow(y-t1[p],2);
       }
std::cout << fehler << std::endl;
   
       return fehler;
}



