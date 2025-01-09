function [geom] = triangolatore(area_ref, marker_triang)

    % global area_ref;
    % global marker_triang;

    if(~exist('mybbtr30.m'))
         addpath('../bbtr30')
         disp('../bbtr30 added to the path')
    end
    
    %----------------------------------------------------------------------------
    %
    % Triangolazione di un dominio quadrata -> triangola un qualsiasi dominio:
    % convesso, concavo, di qualsiasi forma e non è necessario che sia
    % semplicemente connesso, può avere dei buchi
    % con condizioni di Dirichlet sul bordo
    %
    %----------------------------------------------------------------------------
    %
    %  Autore: Stefano Berrone
    %  Politecnico di Torino
    %
    %----------------------------------------------------------------------------
    
    %clc
    %clear all
    
    % -------------------------------
    % Inserimento dei vertici
    % -------------------------------
    
    Domain.InputVertex = [ 0 0
                           1 0
                           1 1
                           0 1];
    
    
    % ---------------------------------------------
    % Definizione del dominio a partire dai Vertici
    % ---------------------------------------------
    
    % Dichiaro le variabili per delimitare il dominio
    Domain.Boundary.Values = 1:4; %Comando che definisce il lati, nel seguente modo:
    % lato di bordo 1 dal nodo 1 al nodo 2
    % lato di bordo 2 dal nodo 2 al nodo 3
    % lato di bordo 3 dal nodo 3 al nodo 4
    % lato di bordo 4 dal nodo 4 al nodo 1
    
    Domain.Holes.Hole = [];       % non ci sono buchi nel dominio
    Domain.Segments.Segment = []; % non ci sono lati forzati nel dominio
    
    % --------------------------------------------------
    % Definizione delle condizioni al contorno a partire
    % dai Vertici e dai lati di bordo
    % --------------------------------------------------
    
    % valori numerici per le condizioni al contorno
    BC.Values = [0.0 12.0 0.0 14.0 0.0 16.0 0.0 0.0 0.0]; % Condizioni di Neumann: 12, 14, 16, 0
    
    % marker delle condizioni al contorno sui bordi del dominio (non il valore che assumono, ma la tipologia), 
    % con la seguente convenzione:
    % dispari -> Dirichlet (si impone sui nodi, il triangolatore aggiunge dei punti sui lati, 
    % che ereditano la condizione di Dirichlet); 
    % pari -> Neumann (si impone sui lati, i punti su quei lati avranno un grado di libertà, 
    % quindi si devono mettere nella matrice)
    
    %BC.Boundary.Values = [3 5 7 9]; % marker dei lati, sono diversi quindi posso imporre condizioni di 
    %Dirichlet diverse sui vari lati
    %BC.Boundary.Values = [2 4 6 8];
    %BC.Boundary.Values = [2 5 7 9];
    BC.Boundary.Values = marker_triang;
    BC.InputVertexValues = [1 1 1 1];% marker dei vertici iniziali, basta che rimangano dispari 
    % Questi indici posso essere anche indici ai valori numerici
    % contenuti nel vettore BC.Values
    
    % Potrei mettere marker dispari e pari sui vari lati a seconda della
    % condizione anche diversa che voglio imporre, i nodi che stanno tra un
    % lato con condizione di Dirichlet e una di Neumann, vince la condizione di
    % Dirichlet
    
    % Il nodo che sta tra due lati di Neumann ha marker 0, come anche tutti i
    % nodi interni, questo vuol dire che hanno 1 grado di libertà
    
    BC.Holes.Hole = [];
    BC.Segments.Segment = [];
    
    
    
    % --------------------------------------------
    % Inserimento dei parametri di triangolazione
    % --------------------------------------------
    
    RefiningOptions.CheckArea  = 'Y'; % controlla che l'area massima sia 0.02 indicata sotto
    RefiningOptions.CheckAngle = 'N'; % check sulla qualità, dovrei specificare l'angolo minimo 
    % ammesso dalla triangolazione, per esempio 20°
    %RefiningOptions.AreaValue  = 0.02; %0.02
    RefiningOptions.AreaValue  = area_ref; 
    RefiningOptions.AngleValue = [];
    RefiningOptions.Subregions = []; % finezze diverse su sottoregioni
    
    
    % --------------------------------------------
    % Creazione della triangolazione e plottaggio
    % --------------------------------------------
    
    [geom] = mybbtr30(Domain,BC,RefiningOptions);
    %draw_grid (geom,1); % plot
    
    % --------------------------------------------------
    % --------------------------------------------------
    
    
    % --------------------------------------------------
    % Rielaborazione dei prodotti del triangolatore
    % per un piu` agevole trattamento delle condizioni
    % al contorno
    % --------------------------------------------------
    
    % Operazioni di pulizia, quando il vettore/matrice creata è più grande del
    % necessario, considero solamente il numero di righe/colonne necessarie
    
    geom.obj.P = geom.obj.P(1:geom.Nobj.N_node,:);
    geom.obj.T = geom.obj.T(1:geom.Nobj.N_ele,:);
    geom.obj.E = geom.obj.E(1:geom.Nobj.N_edge,:);
    geom.obj.Neigh = geom.obj.Neigh(1:geom.Nobj.N_ele,:);
    
    % --------------------------------------------------
    % DOF handler
    
    j  = 1;
    Dj = 1;
    for i=1:size(geom.piv.nlist)
         if geom.piv.nlist(i)==0
            geom.piv.piv(i)=j;
            j = j+1;
         else
            geom.piv.piv(i)=-Dj;
            Dj = Dj + 1;
         end
    end
    
    % --------------------------------------------------
    
    geom.piv.piv = transpose(geom.piv.piv);
    
    % --------------------------------------------------
    
    % geom.pivot.Di dopo le operazioni seguenti contiene l`indice dei nodi
    % di Dirichlet e il corrispondente marker
    
    [X,I] = sort(geom.piv.Di(:,1));
    geom.piv.Di = geom.piv.Di(I,:);
    
    clear X I;
    
    %----------------------------------------------------------------------------

end