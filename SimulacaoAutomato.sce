// Simulador Autômato Estocástico Markoviano com distribuições de fase
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


xf=1;

xs=1;

xt=1;

z= xf+ (xs-1)*h1 + (xt-1)*h2;

NumeroMaximoClientesAtendidos= 10000;  

NClientesTransitorio=400;

NRodadas=10;

Taxa_Abandono=zeros(NRodadas,1)


for  Rodada=1:NRodadas
    
    ClientesAtendidos=0;

    Tempo=0;

    while ClientesAtendidos < NumeroMaximoClientesAtendidos   
           
           TempoAteTransicao=1e10; Evento=-1;
            
            for i = 1:(N+P+2)  // Verifica eventos factíveis
                
                
                if (Af(xf,i)*As(xs,i)*At(xt,i)~=0) then    // evento k é factível em (xf, xs, xt) 
                    
                    Ur= 1e-10+ (1-1e-10)*rand();  // Deslocamento para se eliminar o valor 0;
                    
                    Ti=-(1/Taxas(i))*log(Ur);                    
                    
                    if Ti < TempoAteTransicao then  // Política de corrida (""Race Policy""): escolhe-se o menor tempo
                    
                            TempoAteTransicao= Ti;
                    
                            Evento= i;
                            
                     end// if Ti < TempoAteTransicao then
                   
                    
                end // if (Af(xf,i)*As(xs,i)*At(xt,i)~=0) then 
                             
                    
              end  //  for i = 1:(N+P+2)  // Verifica eventos factíveis  
              
              xf= Af(xf,Evento);
              
              xs= As(xs,Evento);
              
              xt= At(xt,Evento);
              
              z= Af(xf,Evento)+ (As(xs,Evento)-1)*h1 + (At(xt,Evento)-1)*h2;
               
              Tempo=Tempo+ TempoAteTransicao; 
              
               if (Evento==13) then    // Acontecu saída de cliente 
              
                 ClientesAtendidos= ClientesAtendidos +  1;
                 
                  if  (ClientesAtendidos == NClientesTransitorio) then
                   
                            TempoTransitorio= Tempo;
                            
                   end
                 
                 
              end //  if Evento=13 then       
                               
        

      end;  // while

Taxa_Atendimento = (ClientesAtendidos -NClientesTransitorio)/(Tempo-TempoTransitorio)


Taxa_Abandono(Rodada)= Taxas(1)-Taxa_Atendimento 


end  // for  Rodadas=1:NRodadas


// Estimativa de Intervalo de Confiança via Sistribuição t-Student


Estimativa_T_Abandono= mean(Taxa_Abandono);

Estimativa_Var= sum((Taxa_Abandono-Estimativa_T_Abandono).^2 )/(NRodadas-1)

Estimativa_Desvio_Padrao= sqrt(Estimativa_Var/NRodadas);

Tsudent_N_1_95=  2.262;                         // Valor da distribuição t-student para 9 Graus de Liberdade para intervalo 95% de confiança

Delta= Tsudent_N_1_95* Estimativa_Desvio_Padrao;   

Intervalo=zeros(2,1);

Intervalo(1)= Estimativa_T_Abandono -Delta;
Intervalo(2)= Estimativa_T_Abandono + Delta;


DuracaoSimulacao=toc();  // Término cronometro


        

