function bandplot3D(elementos,nodos,variable,lims,lineColor,nColores,minimalist, faceAlpha, lineStyle, colorbarTitle)

% BANDPLOT  Graficador de variables.
% 
% BANDPLOT(elementos,nodos,variable)
% BANDPLOT(elementos,nodos,variable,lims)
% BANDPLOT(elementos,nodos,variable,lims,lineColor)
% BANDPLOT(elementos,nodos,variable,lims,lineColor,nColores)
% 
% Numeración de los nodos de los elementos:
%  4---7---3
%  |       |
%  8   9   6
%  |       |
%  1---5---2
% 
% elementos: Matriz de conectividades.
% nodos:     Matriz de coordenadas nodales.
% variable:  Matriz m x n con la variable a graficar, donde [m n] = size(elementos).
% lims:      Vector [cmin cmax] con los límites a utilizar para la variable.
% lineColor: String que especifica el color para los bordes de los
%            lementos. Utilizar 'none' para no mostrar los bordes.
% nColores:  Cantidad de bandas de colores a utilizar.
%minimalist: No node nor element numbering
Nodperelement = size(nodos,2);

% error(nargchk(3, 6, nargin))
isValidVariable = all(size(variable) == size(elementos));
if ~isValidVariable
    error('size(elementos) debe ser igual a size(variable)')
end
    
if (nargin < 4) || isempty(lims)
    lims = [min(min(variable))-0.001 max(max(variable))+0.001];
end

if (nargin < 5) || isempty(lineColor)
    lineColor = 'none';
end

tol = 1E-5;
nMinCol = 3;
if (nargin < 6) || isempty(nColores)
    if diff(lims) < lims(2)*tol
        nColores = nMinCol;
    else
        nColores = 10;
    end
end

if (nargin < 7) || isempty(minimalist)
    minimalist = false;
end

if (nargin < 8) || isempty(faceAlpha)
    faceAlpha = 0.5;
end

if (nargin < 9) || isempty(lineStyle)
    lineStyle = "-";
end

if (nargin < 10) || isempty(colorbarTitle)
    colorbarTitle = '';
end

if diff(lims) < lims(2)*tol
    lims = lims + [-1 1]*lims(2)*tol;
end

if nColores < nMinCol
    nColores = nMinCol;
end

nNodos = size(elementos,2);
nel = size(elementos,1);

switch nNodos
    case {8}
        vNod1 = [1 2 3 4];
        vNod2 = [1 2 6 5];
        vNod3 = [2 6 7 3];
        vNod4 = [1 5 8 4];
        vNod5 = [5 6 7 8];
        vNod6 = [4 3 7 8];
    case {20}
       vNod1 = [1 9 2 10 3 11 4 12];
%        vNod2 = [2 18 6 14 7 19 3];
       vNod2 = [2 18 6 14 7 19 3 10];
       vNod3 = [4 12 1 17 5 16 8 20];
       vNod4 = [5 13 6 14 7 15 8 16];
       vNod5 = [1 9 2 18 6 13 5 17];
       vNod6 = [4 11 3 19 7 15 8 20];
%        vNod2 = [2 14 6 18 7 15 3];
%        vNod3 = [4 12 1 13 5 20 8 16];
%        vNod4 = [5 17 6 18 7 19 8 20];
%        vNod5 = [1 9 2 14 6 17 5 13];
%        vNod6 = [4 11 3 15 7 19 8 16];
    otherwise
        vNod = 1:nNodos;
end

% Creo los patches
for iele = 1:nel
    eleNodes = elementos(iele,:);
    h = patch('Faces',vNod1,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    h = patch('Faces',vNod2,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    h = patch('Faces',vNod3,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    h = patch('Faces',vNod4,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    h = patch('Faces',vNod5,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    h = patch('Faces',vNod6,'Vertices',nodos(eleNodes,:),'FaceVertexCData',variable(iele,:)');
    set(h,'FaceColor','interp','EdgeColor',lineColor,'lineStyle',lineStyle,'CDataMapping','scaled','FaceAlpha',faceAlpha);
    
    xx = nodos(eleNodes,1);
    yy = nodos(eleNodes,2);
    zz = nodos(eleNodes,3);
    
    if minimalist == false
        text(mean(xx),mean(yy),mean(zz),num2str(iele),'VerticalAlignment','bottom','Color','b','FontSize',8);
    end
end
for i = 1:size(nodos,1)
    xx(i) = nodos(i,1);
    yy(i) = nodos(i,2);
    zz(i) = nodos(i,3);
    
    if minimalist==false
        text(xx(i),yy(i),zz(i),num2str(i),'VerticalAlignment','bottom','Color','r','FontSize',8);
    end
end
% Acomodo escalas y demás
colormap(jet(nColores))
caxis(lims)

if nColores <= 20
    nTicks = nColores;
else
    nTicks = 20;
end

ticks = lims(1):((diff(lims))/nTicks):lims(2);
tickLabels = cell(size(ticks));

for iTick = 1:length(ticks)
    tickLabels{iTick} = sprintf('%6.5E',ticks(iTick));
end
ax = gca;
hcb = colorbar('YTick',ticks,'YTickLabel',tickLabels);
colorTitleHandle = get(hcb,'Title');
set(colorTitleHandle ,'String',colorbarTitle,'Interpreter','latex', 'fontSize',14);

view(3)
set(ax,'XTick',[],'YTick',[],'ZTick',[],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'visible','on')
ax.XLabel.String = "X";
ax.YLabel.String = "Y";
ax.ZLabel.String = "Z";
daspect([1 1 1])
box on
axis auto

