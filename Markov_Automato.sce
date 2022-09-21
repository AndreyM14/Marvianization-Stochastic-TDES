// Conversor Autômato Estocástico Non Markoviano em Cadeia de Markov: Casamendo da média e da variância 
// Prof. Carlos Andrey Maia
// Departamento de Eng. Elétrica da UFMG
//       30/11/2021

// ---------------   Entradas  ----------------------

// X_i={ 1, 2, ...,gi}

// E= E1 U E2  U   ... En = { 1, ..., L}  união de todos conjunto dos eventos

// A_i: Matriz de incidencia  do autômato Gi x_i^+ =  A_i(xi,ei) =fi(xi,ei), xi \in X_i e ei \in E

// A_i(xi,ei) = 0 se ei não é factível no estado xi;  A_i(xi,ei) = xi se ei não pertence a Ei.


// Simulação de Fila com atendimento não markoviano com abandono de clientes

// Conjunto de evento  considerados    E = { a, c, q1, q2, ..., qn, d, f1, ...f_{P-1}}

//  a : chegada de cleientes à fila; 
// c: admissão de cliente pelo servidor;
// qi: indica que clinetes abandonam  da fila qdo o comprimento da fila é i
// d: indica que cliente foi atendido com êxito pelo servidor
// fi: eventos intermediários que indicam saltos markovianos ( são utilizados para representação de tempo de serviços não marckovianos)


// Estados para o espaço de espera: Xf = [ 1, 2, ..., N,N+1}   ] 


//  Estados para o servidor: Xs = [ 1, 2}   ]  


// Estados para a temporização:  Xf = [ 1, 2, ...,P}   ]   

//

clear;

tic();  // Inicio cronometro


N=10;  P=4;;

Af= ones(N+1, N+P+2 ) ;

for j=2:(N+1)
    Af(j,:)=j* Af(j,:);  // Inicia-se com auto-laços para todos os eventos
end

Af(:, 1: N+2)= 0*Af(:, 1: N+2); // Preenche com zeros eventos do alfabeto

Af(:,1)= [ 2 3 4 5  6 7 8 9 10 11 0 ]';  // Evento a: chegada de clientes

Af(:, 2) = [ 0  1  2  3  4  5  6  7  8  9  10];  // Evento c: admissão de clientes pelo servidor

Af(2,3)=1;  // Evento q1 : abandono da fila quando há um cliente esperando
Af(3,4)=2 ;  // Evento q2  : " ....                   dois clientes ...."
Af(4,5)=3 ;  // ...
Af(5,6)= 4;  
Af(6,7)= 5;
Af(7,8) =6;
Af(8,9)=7;   // ...
Af(9,10)=8;
Af(10,11)=9;
Af(11,12)=10;  // Evento qn:  bandono da fila quando há n clientes esperando 


// --------------- Autômato para o servidor  ---------------------------

// XS = [ 1, 2}   ]  E = { a, c, q1, q2, ..., qn, d, f1, ...f_{P-1}}


As= ones(2, N+P+2 );

As(2,:)= 2* As(2,:);

As(:,2)=[2 0]';  //  Evento c

As(:,2+N+1)=[0 1]';  //  Evento d

As(:,2+N+2)=[0 2]';  //  Evento f1


// Auômato com  saltos markovianos para o tempo de serviço

// XS = [ 1, 2, ...,P}   ]  E = { a, c, q1, q2, ..., qn, d, f1, ...f_{P-1}}

At= ones(P, N+P+2 );

At(2,:) =2*At(2,:); At(3,:) =3*At(3,:);  At(4,:) =4*At(4,:);   // Auto-laços

At(:, N+3:N+6)=0*At(:, N+3:N+6);  // Zero indica evento não factivel

At(1,N+4)=2;  // evento f1

At(2,N+5)=3;  // evento f2

At(3,N+6)=4;  // evento f3

At(4, N+3)= 1;  // evento d


//------------Codificaçção vetorial dos estados

// h1= N+1;  h2= h1*2

// (xf, xs, xt)  <--->  k= xf+ (xs-1)*h1 + (xt-1)*h2


h1= N+1;  h2= h1*2;


// -----------------------------------Taxas para o saltos markovianos  ---------

Taxas=zeros(1, (N+P+2));

Taxas(1)= 2; // Evento a   //2

Taxas(2)= 10;  // Evento c  //10


for i=1:10
    
    Taxas(i+2)= i*0.2;  // Evento qi
    
end

Taxas(13)=4*2.5;  // Evento d

Taxas(14) =  4*2.5;  // Evento f1

Taxas(15)=  4*2.5;  //  Evento f2

Taxas(16)=  4*2.5;  // Evento f3


//---------------Inicializaçãode matrizes--------------------------------------

Nx= (N+1)*2*P;

Input=zeros(Nx,(N+P+2));

Output=zeros(Nx,(N+P+2));

Aii=zeros(Nx,1);

Am=zeros(Nx,Nx);

Acessivel=zeros(Nx,1);    Acessivel(1)=1;  // Estados acessíveis



// -----------------------------------------------------------------------------

for xt = 1:P   // Estados para as fases
 
    for xs = 1:2  // Estados do servidor 
        
        for xf = 1:(N+1)  // Estados do espaçoe de espera 
                  
            z= xf+ (xs-1)*h1 + (xt-1)*h2
            
            for i = 1:(N+P+2)  // Verifica eventos factíveis
                
                if (Af(xf,i)*As(xs,i)*At(xt,i)~=0) then    // evento k é factível em (xf, xs, xt) 
                    
                    zi= Af(xf,i)+ (As(xs,i)-1)*h1 + (At(xt,i)-1)*h2   // z é um estado de entrada para zi para o evento i
                    
                    Input(zi,i)= z;  // f(z,i) = zi
                    
                    Output(z,i)=zi;
                    
                    // AcessZ=  Acessivel(z)  
                    
                    //AcessZi= max( Acessivel(zi),Acessivel(z) )  
                    
                    Acessivel(zi)=  max( Acessivel(zi),Acessivel(z) )                  
                    
                    Am(zi,z)=Taxas(i);  // Matriz completa
                    
                    Aii(z)= Aii(z) - Taxas(i);    // Aii(zi)*zi  + Somatorio (Input(zi,i)*Taxas(i))=0   // Equilíbrio do Sitemas Markoviano
                                      
                end // if 
                
                Am(z,z)= Aii(z);
                    
              end                      
        
         end;
        
    end;
end



// -------------  Estados alcançáveis   : Modelo Tropical


Alcancavel=-1e100*ones(Nx,1); Alcancavel(1)=1;


for k=1:Nx    
    
  for z=1:Nx
    
    for i = 1:(N+P+2)
        
        zin= Input(z,i);
        
        if    zin >0 then
            
        Alcancavel(z)= max(Alcancavel(z), Alcancavel(zin) )
         
        end // if
           
            
    end // for
    
  end // for  z=1:Nx

end  // 1:Nx  
    

// -------------------   Condição inicial para o Algorítmo

N_nulo=0

xi=ones(Nx,1);

for k=1:Nx 
    
    if Alcancavel(k) < 1   then
        
        N_nulo= N_nulo +1;
        
        xi(k)=0;
        
     end // if
                  
end  

xi= xi/(Nx-N_nulo);





// Solução do sistema Amodif: Amodif*x=b,   Amodif = [ 1  1 1 ....1; A(2:n,:)]  b = [ 1  0 0  .... 0]']  ;  Método de Gauss-Seidel


//----------------------------------------------------Metódo completo 1----------------------------
//b=-Am(2:Nx, 1)  // Fixando 1 como solução

//y=inv(Am(2:Nx, 2:Nx))*b

//p=[ 1; y]  // Solução não normalizada

// Pn=p/sum(p)

//     -------------------------------------- Solução pela  inversão da Matriz  ---------------------------------------


//Am2=[ones(1,Nx); Am(2:Nx,:)];  // Remove-se a primeira linha da matriz

//b=zeros(Nx,1); b(1)=1;

//Pn= inv(Am2)*b;

// Erro_metodo_completo= sum(abs(Am*Pn));


// ----------------- Solução pelo Método de Jacobi com Filtragem de Oscilações ------------------------------


// x=zeros (Nx, 1); x(1)=1;
 
//  x=Pn;
 
x=xi

ErroMax=10;

x_hat=xi;

alpha=0.9;

while abs(ErroMax) > 0.0001
    
    ErroMax=0;
    
    xa=x;
    
    x(1) = (1 - sum(xa(2:Nx)) );
    
       
    ErroMax= max(ErroMax,  abs(x(1)-xa(1))/xa(1))
    
    
    for z=2:Nx
        
        Soma_Lambdai_x_i=0;
                  
        if Alcancavel(z) > 0.5   then
        
            for i = 1:(N+P+2)
            
                  zi= Input(z,i);  // z=f(zi,i) 
           
                  if zi >0 then
               
                      Soma_Lambdai_x_i = Soma_Lambdai_x_i+ Taxas(i)*xa(zi);  //  
            
                  end // if
                              
                x_hat(z) = (- Soma_Lambdai_x_i )/ Aii(z);
                
                x(z)= alpha* x_hat(z)+(1-alpha)*xa(z)  // Filtro  para  oscilações muito grandes
     
             end  // For
        
        end // if if Alcancavel(k) > 0 
            
       ErroMax= max(ErroMax, abs(x(z)-xa(z))/xa(z) );
          
   end  //  for z=2:Nx
 
 
end  // While

// MaximoErroMetodo= max(abs(Pn-x)./Pn)



// Taxa efeitva de clientes atendidos

 Lambda_ef=0;

 for xf = 1:(N+1),
     
      z= xf+ (2-1)*h1 + (4-1)*h2
      
      Lambda_ef= Lambda_ef+ x(z)*Taxas(13)
      
end


PercentualAbandono = (Taxas(1)- Lambda_ef)/Taxas(1)


DuracaoSimulacao=toc();  // Término cronometro
  
  







