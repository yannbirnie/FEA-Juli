clear; clc; close all

escala=1; %Escala desplazamientos

CASO = "corte" % "corte", "traction" 
%% Discretizacion

%H8
x_pos = [ 0 ,(1-sqrt(0.6))/2 , 1-(1-sqrt(0.6))/2 , 1];
y_pos = x_pos;
z_pos = [0 1/3 2/3 1];

nodes = nan(4^3,3);
counter = 0;
for z = z_pos
    for y = y_pos
        for x = x_pos
            counter = counter+1;
            nodes(counter,:) = [x y z];
        end
    end
end

elements = nan(27,8);
element_1 = [1 5 6 2 17 21 22 18];
counter = 0;
for kk = 0:2 % Creating the neighbouring element following the z axis
    for jj = 0:2 % Creating the neighbouring element following the y axis
        for ii = 0:2 % Creating the neighbouring element following the x axis
            counter=counter+1;
            elements(counter,:)= element_1  + ii + jj*4 + kk*16 ;
        end
    end
end

%% Propiedades del Material
E = 5; %E=207e3; %CHECK UNIDADES [MPa]
nu=0.3;
G=E/2/(1+nu);

%relacion constitutiva
C = (E/((1+nu)*(1-2*nu)))*[1-nu    nu   nu      0   0   0
                           nu    1-nu   nu      0   0   0
                           nu      nu  1-nu     0   0   0
                            0   0   0   0.5-nu      0   0
                            0   0   0   0   0.5-nu  0
                            0   0   0   0   0   0.5-nu];
                        
%% Definiciones
nDofNod = 3;                    % grados de libertad por nodo
nel = size(elements,1);         % elementos
nNod = size(nodes,1);           % nodos
nNodEle = size(elements,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad
nDims = size(nodes,2);          % dimensiones del problema

 bandplot3D(elements,nodes, zeros(size(elements)),[-100 100], 'k');

%% Puntos de Gauss
[wpg, upg, npg] = gauss([3 3 3]);
%generaciôn de funciones de forma y sus derivadas evaluadas en los puntos de gauss
N = zeros(1,nNodEle,npg);    dN = zeros(3,nNodEle,npg);
for ipg = 1:npg
     % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        x = ksi;    y = eta;    z = zeta;
        
        % Funciones de forma respecto de ksi, eta,zeta del H8 !!!
        N(:,:,ipg) = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];

        % Derivadas de las funciones de forma respecto de ksi, eta & zeta
        % para el H8
        dN(:,:,ipg) = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                        x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                        x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];
end

%% Matriz de rigidez
K = zeros(nDofTot,nDofTot);
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN(:,:,ipg)*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,ipg);          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:nDofNod*nNodEle-2)  = dNxy(1,:);
        B(2,2:3:nDofNod*nNodEle-1)  = dNxy(2,:);
        B(3,3:3:nDofNod*nNodEle)    = dNxy(3,:);
        B(4,1:3:nDofNod*nNodEle-2)  = dNxy(2,:);
        B(4,2:3:nDofNod*nNodEle-1)  = dNxy(1,:);
        B(5,2:3:nDofNod*nNodEle-1)  = dNxy(3,:);
        B(5,3:3:nDofNod*nNodEle)    = dNxy(2,:);
        B(6,1:3:nDofNod*nNodEle-2)  = dNxy(3,:);
        B(6,3:3:nDofNod*nNodEle)    = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
    end
    eleDofs = node2dof(elements(iele,:),nDofNod); 
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;
end

%% Condiciones de borde
bc = false(nNod,nDofNod);       % Matriz de condiciones de borde

bc(1:16,[3]) = true;
bc([1],[1 2 3]) = true;
bc([13],[1]) = true;
    
%%
% figure(1)
% MeshplotTrigMec(elementos,nodos,bc,'k',1)

%% Puntos de Gauss superficiales
[wpg, upg, npg] = gauss([3 3]); 

%% Cargas de  presión hydrostática en Z
presion = 1; %MPa ---> Presion aplicada sobre la cara

R = zeros(nNod,nDofNod);  
area = 0;
if CASO == "traction"
    for iele = [19 20 21 22 23 24 25 26 27]
        areaelemental = 0;
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss superficiales
            ksi = upg(ipg,1);
            eta = upg(ipg,2);
            zeta = 1;
           
            x = ksi;    y = eta; z = zeta;

            % Funciones de forma respecto de ksi, eta,zeta del H8 !!!
            N = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];
            dN = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                   x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                   x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

            jac = dN*nodesEle;
            areaelemental = areaelemental + wpg(ipg)*cross(jac(2,:),jac(3,:));

            Re = Re + N'*wpg(ipg)*cross(jac(1,:),jac(2,:))*presion; %aplicando 1Mpa de presiôn sobre la cara trasera en direcciôn x

        end

        area = area + areaelemental;
        R(elements(iele,:),:) = R(elements(iele,:),:) + Re;
    end
end 
%% Cargas corte en x
if CASO == "corte"
    R = zeros(nNod,nDofNod);     
    for iele = [1 2 3 4 5 6 7 8 9]
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss superficiales

                ksi = upg(ipg,1);
                eta = upg(ipg,2);
                zeta = -1;
                x = ksi;    y = eta;    z = zeta;

                carga = [1 0 0
                         1 0 0
                         1 0 0
                         1 0 0
                         0 0 0
                         0 0 0
                         0 0 0
                         0 0 0];

                %  Funciones de forma respecto de ksi, eta,zeta del H8 !!!
                N = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];
                ENE = [N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0 0
                       0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0
                       0 0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8)];

    %             ENE = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)];

                dN = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                       x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                       x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

                jac = dN*nodesEle;

                Re = Re + reshape(ENE'*ENE*reshape(carga',[],1)*norm(cross(jac(1,:),jac(2,:)))'*wpg(ipg),3,[])'; %aplicando 1Mpa de presiôn sobre la cara superior en direcciôn x
         end
        R(elements(iele,:),:) = R(elements(iele,:),:) + Re;
    end
    for iele = [19 20 21 22 23 24 25 26 27]
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss superficiales

                ksi = upg(ipg,1);
                eta = upg(ipg,2);
                zeta = 1;
                x = ksi;    y = eta;    z = zeta;

                carga = [0 0 0
                         0 0 0
                         0 0 0
                         0 0 0
                         -1 0 0
                         -1 0 0
                         -1 0 0
                         -1 0 0];

             %  Funciones de forma respecto de ksi, eta,zeta del H8 !!!
                N = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];
                ENE = [N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0 0
                       0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0
                       0 0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8)];
    %             ENE = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)];
                dN = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                       x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                       x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

                jac = dN*nodesEle;

               Re = Re + reshape(ENE'*ENE*reshape(carga',[],1)*norm(cross(jac(1,:),jac(2,:)))'*wpg(ipg),3,[])'; %aplicando 1Mpa de presiôn sobre la cara inferior en direcciôn -x
         end
        R(elements(iele,:),:) = R(elements(iele,:),:) + Re;
    end
    for iele = [3 6 9 12 15 18 21 24 27]
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss superficiales
                ksi = 1;
                eta = upg(ipg,1);
                zeta = upg(ipg,2);
                x = ksi;    y = eta;    z = zeta;

                carga = [0 0 0
                         0 0 0
                         0 0 -1
                         0 0 -1
                         0 0 0
                         0 0 0
                         0 0 -1
                         0 0 -1];

             %  Funciones de forma respecto de ksi, eta,zeta del H8 !!!
                N = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];
                ENE = [N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0 0
                       0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0
                       0 0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8)];
    %             ENE = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)];
                dN = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                       x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                       x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

                jac = dN*nodesEle;

               Re = Re + reshape(ENE'*ENE*reshape(carga',[],1)*norm(cross(jac(2,:),jac(3,:)))'*wpg(ipg),3,[])'; %aplicando 1Mpa de presiôn sobre la cara ksi=1 en direcciôn -z
         end
        R(elements(iele,:),:) = R(elements(iele,:),:) + Re;
    end
    for iele = [1 4 7 10 13 16 19 22 25]
        Re = zeros(nNodEle,nDofNod);
        nodesEle = nodes(elements(iele,:),:);
        for ipg = 1:npg
            % Punto de Gauss superficiales
                ksi = -1;
                eta = upg(ipg,1);
                zeta = upg(ipg,2);
                x = ksi;    y = eta;    z = zeta;

                carga = [0 0 1
                         0 0 1
                         0 0 0
                         0 0 0
                         0 0 1
                         0 0 1
                         0 0 0
                         0 0 0];

             %  Funciones de forma respecto de ksi, eta,zeta del H8 !!!
                N = [ (x*y)/8 - y/8 - z/8 - x/8 + (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, y/8 - x/8 - z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, x/8 + y/8 - z/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8, x/8 - y/8 - z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, z/8 - y/8 - x/8 + (x*y)/8 - (x*z)/8 - (y*z)/8 + (x*y*z)/8 + 1/8, y/8 - x/8 + z/8 - (x*y)/8 - (x*z)/8 + (y*z)/8 - (x*y*z)/8 + 1/8, x/8 + y/8 + z/8 + (x*y)/8 + (x*z)/8 + (y*z)/8 + (x*y*z)/8 + 1/8, x/8 - y/8 + z/8 - (x*y)/8 + (x*z)/8 - (y*z)/8 - (x*y*z)/8 + 1/8];
                ENE = [N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0 0
                       0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8) 0
                       0 0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0 N(5) 0 0 N(6) 0 0 N(7) 0 0 N(8)];
    %             ENE = [N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)
    %                    N(1) N(2) N(3) N(4) N(5) N(6) N(7) N(8)];
                dN = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                       x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                       x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

                jac = dN*nodesEle;

               Re = Re + reshape(ENE'*ENE*reshape(carga',[],1)*norm(cross(jac(2,:),jac(3,:)))'*wpg(ipg),3,[])'; %aplicando 1Mpa de presiôn sobre la cara ksi=-1 en direcciôn +z
         end
        R(elements(iele,:),:) = R(elements(iele,:),:) + Re;
    end
end
%% Plot de las cargas
hold on 
scale = 2;
quiver3(nodes(:,1),nodes(:,2),nodes(:,3),R(:,1),zeros(size(R(:,2))),zeros(size(R(:,3))), scale ,'b');
quiver3(nodes(:,1),nodes(:,2),nodes(:,3),zeros(size(R(:,1))),R(:,2),zeros(size(R(:,3))), scale ,'y');
quiver3(nodes(:,1),nodes(:,2),nodes(:,3),zeros(size(R(:,1))),zeros(size(R(:,2))),R(:,3), scale ,'r');

 
 %% ordenamiento vectorial
isFixed = reshape(bc',[],1);
isFree = ~isFixed;
Rr = reshape(R',[],1);

% Solver
D = zeros(nDofTot,1);
D(isFree) = K(isFree,isFree)\Rr(isFree);
Rr(isFixed) = K(isFixed,isFree)*D(isFree);

% Reacciones

reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rr(isFixed);
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Recuperación de tensiones en los nodos
dN = zeros(3,nNodEle, nNodEle);
%Local node position of the H8
uNod = [-1  -1  -1
        -1   1  -1
         1   1  -1
         1  -1  -1
        -1  -1   1
        -1   1   1
         1   1   1
         1  -1   1];

     %calculo de las derivadas de las funciones de forma en los unods
     for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        zeta = uNod(inode,3);
        x=ksi;  y=eta;  z=zeta;
  %     Derivadas de las funciones de forma respecto de ksi, eta & zeta
  %     para el H8
        dN(:,:,inode) = [ y/8 + z/8 - (y*z)/8 - 1/8, z/8 - y/8 + (y*z)/8 - 1/8, y/8 - z/8 - (y*z)/8 + 1/8, (y*z)/8 - z/8 - y/8 + 1/8, y/8 - z/8 + (y*z)/8 - 1/8, - y/8 - z/8 - (y*z)/8 - 1/8, y/8 + z/8 + (y*z)/8 + 1/8, z/8 - y/8 - (y*z)/8 + 1/8
                          x/8 + z/8 - (x*z)/8 - 1/8, (x*z)/8 - z/8 - x/8 + 1/8, x/8 - z/8 - (x*z)/8 + 1/8, z/8 - x/8 + (x*z)/8 - 1/8, x/8 - z/8 + (x*z)/8 - 1/8, z/8 - x/8 - (x*z)/8 + 1/8, x/8 + z/8 + (x*z)/8 + 1/8, - x/8 - z/8 - (x*z)/8 - 1/8
                          x/8 + y/8 - (x*y)/8 - 1/8, x/8 - y/8 + (x*y)/8 - 1/8, - x/8 - y/8 - (x*y)/8 - 1/8, y/8 - x/8 + (x*y)/8 - 1/8, (x*y)/8 - y/8 - x/8 + 1/8, y/8 - x/8 - (x*y)/8 + 1/8, x/8 + y/8 + (x*y)/8 + 1/8, x/8 - y/8 - (x*y)/8 + 1/8];

     end
stress = zeros(nel,nNodEle,6);
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Derivadas de x,y, respecto d % '2'; % e ksi, eta
        jac  = dN(:,:,inode)*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN(:,:,inode);          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle); %dimensions 6x(3*20)
        B(1,1:3:nDofNod*nNodEle-2)  = dNxy(1,:);
        B(2,2:3:nDofNod*nNodEle-1)  = dNxy(2,:);
        B(3,3:3:nDofNod*nNodEle)    = dNxy(3,:);
        B(4,1:3:nDofNod*nNodEle-2)  = dNxy(2,:);
        B(4,2:3:nDofNod*nNodEle-1)  = dNxy(1,:);
        B(5,2:3:nDofNod*nNodEle-1)  = dNxy(3,:);
        B(5,3:3:nDofNod*nNodEle)    = dNxy(2,:);
        B(6,1:3:nDofNod*nNodEle-2)  = dNxy(3,:);
        B(6,3:3:nDofNod*nNodEle)    = dNxy(1,:);
        
        eleDofs = node2dof(elements(iele,:),nDofNod);
%       stress(inode,:,iele) = C*B*D(eleDofs);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end
    
%% Configuración deformada
figure
nodePosition = nodes + escala*(reshape(D,nDofNod,[]))';

bandplot3D(elements,nodes, zeros(size(elements)),[-100 100], 'k', 3, true, 0.1, '--');
bandplot3D(elements,nodePosition, zeros(size(elements)),[-100 100], 'k', 3, true);

% figure(2)
% % meshplot(elements,nodePosition,bc,'r',0)
% meshplot(elements,nodePosition,'r',0)

%% Gráfico
figure
% bandplot3D(elements,nodePosition,stress(:,:,1),[],'k');
% bandplot3D(elements,nodePosition,stress(:,:,1),[],'k');
if CASO == "traction"
    bandplot3D(elements,nodes+ escala*(reshape(D,nDofNod,[]))',stress(:,:,3),[],'k',3,true, 0.5, '-', '$\sigma_{zz}$ [MPa]');
elseif CASO == "corte"
    bandplot3D(elements,nodes+ escala*(reshape(D,nDofNod,[]))',stress(:,:,6),[],'k',1,true, 0.5, '-', '$\sigma_{xz}$ [MPa]');
end