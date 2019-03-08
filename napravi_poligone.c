#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    
    int x0,y0; //centar
    int r;
    int n;
    
    fprintf(stderr,"Unesi x0, y0, r i broj tacaka n poligona P\n");
    scanf("%d %d %d %d",&x0,&y0,&r,&n);
    
    printf("%d\n",n);
    for(int j=0;j<n;j++){
        int x = (int)(x0+r*cos((j*2*M_PI)/n));
        int y = (int)(y0+r*sin((j*2*M_PI)/n));
        
        printf("%d %d\n",x,y);
        
    }
    //printf("%d %d\n",(int)(x0+r*cos((0*2*M_PI)/n)),(int)(y0+r*sin((0*2*M_PI)/n)));
    
    //printf("\n");
    
    fprintf(stderr,"Unesi x0, y0, r i broj tacaka m poligona Q\n");
    scanf("%d %d %d %d",&x0,&y0,&r,&n);
    
    printf("%d\n",n);
    for(int j=0;j<n;j++){
        int x = (int)(x0+r*cos(1.27+(j*2*M_PI)/n));
        int y = (int)(y0+r*sin(1.27+(j*2*M_PI)/n));
        
        printf("%d %d\n",x,y);
        
    }
    
    //printf("%d %d\n",(int)(x0+r*cos(1.27+(0*2*M_PI)/n)),(int)(y0+r*sin(1.27+(0*2*M_PI)/n)));
       
    
    exit(EXIT_SUCCESS);
}
