function [geom] = triangolatore(area_ref, marker_triang, marker_vertici)


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
    BC.InputVertexValues = marker_vertici;% marker dei vertici iniziali, basta che rimangano dispari 
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

    %% Parte 1: Aggiungo le informazioni geometriche alla struttura
    nnode = geom.Nobj.N_node;
    
    for l=1:geom.Nobj.N_edge
    
      n(1)=geom.obj.E(l,1); % edge starting node
      n(2)=geom.obj.E(l,2); % edge ending node
      e(1)=geom.obj.E(l,3); % side element
      e(2)=geom.obj.E(l,4); % side element
      
      nnode = nnode + 1;
    
      geom.obj.P(nnode,:) = (geom.obj.P(n(1),:)+...
					    geom.obj.P(n(2),:))/2; % P of the edge midpoint
    
      geom.obj.E(l,5)=nnode; % to connect the edge with its midpoint
    
      % Ricerca della posizione del punto medio
      idx = [1 2 3];
    
      for el=e
        
        if(el ~= -1)
          acc = 0;
          acc = idx * ( geom.obj.T(el,1:3)==n(1) )'; % Ricerco la posizione dell'indice globale e la converto in indice 
          % locale, moltiplicando scalarmente per idx
          acc = acc + idx * ( geom.obj.T(el,1:3)==n(2) )'; % Il punto medio sta tra:
          % 1 e 2 -> somma 3 -> posizione 4
          % 1 e 3 -> somma 4 -> posizione 5
          % 2 e 3 -> somma 5 -> posizione 6
    
          switch acc
	    case 3
	      geom.obj.T(el,4) = nnode;
	    case 4
	      geom.obj.T(el,6) = nnode;
	    case 5
	      geom.obj.T(el,5) = nnode;
	    otherwise
	      disp('sconoscuto');
          end % switch acc      
        end
    
      end % for el=e
    
    %% Parte 2: Distinguo se il lato di bordo ha la condizione di Neumann o di Dirichlet
    % Se ho un lato che non è stato tagliato dal triangolatore e i due vertici
    % hanno D, non è detto che il punto medio sia D, se su quel lato ho N devo
    % imporre N anche nel punto medio-> devo vedere se quello è un lato
    % iniziale e se ha la condizione di N, in questo modo impongo la condizione
    % giusta al punto medio
    % Obiettivo: non voglio fare ricerche all'interno di vettori per es. -> log(n)
    
      VertexValue = [0 0];
    %  Vertex = [0 0];
      D = [0 0];
      InputVertexValue=[0 0];
      
    %idxV = 1:length(geom.input.BC.InputVertexValues); %indice dei vertici
    
				    % Se il lato e` di bordo
      if( any(e==-1) )
				    % Lato di Dirichlet
        
        if( geom.piv.nlist(n(1))~=0 && geom.piv.nlist(n(2))~=0 )
				    %-----------------------
          if( geom.piv.nlist(n(1)) ~= geom.piv.nlist(n(2)) )
				    %
	    if( any(geom.input.BC.InputVertexValues==geom.piv.nlist(n(1))) )
				    % Vertice(1) con marker speciale
	      VertexValue(1) = 1;
		    % e` il vettore con un 1 in corrispondenza del vertice
		    % del lato con marker speciale
	      
    %	  Vertex(1) = geom.input.BC.Boundary.Values*...
    %		      (geom.input.BC.InputVertexValues == geom.piv.nlist(n(1)))';
                 % valore del marker del lato del poligono iniziale che segue n1
	      
	         InputVertexValue(1) = [1:length(geom.input.BC.InputVertexValues)]*...
			         (geom.input.BC.InputVertexValues == geom.piv.nlist(n(1)))';
                 % marker del nodo del poligono iniziale che corrisponde al nodo n(1) del mio lato
    
	          % valore del marker del vertice del poligono iniziale n2
	      D(1) = geom.piv.nlist(n(2));
	    end
				    %	
	    if( any(geom.input.BC.InputVertexValues==geom.piv.nlist(n(2))) )
	      VertexValue(2) = 1;
    %	  Vertex(2) = geom.input.BC.Boundary.Values*...
    %		      (geom.input.BC.InputVertexValues == geom.piv.nlist(n(2)))';
              InputVertexValue(2) = [1:length(geom.input.BC.InputVertexValues)]*...
			         (geom.input.BC.InputVertexValues == geom.piv.nlist(n(2)))';
	      D(2) = geom.piv.nlist(n(1));
	    end
				    %
	    if( sum(VertexValue) ~= 2 )
	      Di = VertexValue*D';
				    % nodo con condizione di Dirichlet
              geom.piv.nlist(nnode)= Di;
              geom.piv.Di(end+1,:) = [nnode, Di];
              geom.piv.piv(nnode) = min(geom.piv.piv)-1;
	    else
                  % diamo al nuovo nodo il marker del lato
                  % il lato che stiamo analizzando e` un lato del poligono
                  % iniziale:
              
	     % l'indice del lato e` quello del nodo di inizio di quel lato
	      if( max(InputVertexValue)-min(InputVertexValue)>1 ) % siamo sul lato di chiusura
	        Di = geom.input.BC.Boundary.Values(max(InputVertexValue));
	      else % siamo sui lati 1->2->3->4->
	        Di = geom.input.BC.Boundary.Values(min(InputVertexValue));
	      end
			        % check della condizione di Neumann aperta
	      if( rem(Di,2)== 0 ) % nodo con grado di liberta`, lato di
                                  % Dirichlet aperto
	        geom.piv.nlist(nnode)= 0;
	        geom.piv.piv(nnode) = max(geom.piv.piv)+1;
                disp('non dovevi essere qui');
	      else
                geom.piv.nlist(nnode)= Di;
                geom.piv.Di(end+1,:) = [nnode, Di];
                geom.piv.piv(nnode) = min(geom.piv.piv)-1;
	      end
	    end % if( sum(VertexValue) ~= 2 )
	        %----------------------------------
          else % if( geom.piv.nlist(n(1)) ~= geom.piv.nlist(n(2)) )
	    Di = geom.piv.nlist(n(1));
            geom.piv.nlist(nnode)= Di;
            geom.piv.Di(end+1,:) = [nnode, Di];
            geom.piv.piv(nnode) = min(geom.piv.piv)-1;
          end % if( geom.piv.nlist(n(1)) ~= geom.piv.nlist(n(2)) )
	      %----------------------------------
    
        else
				    % Lato di Neumann
          geom.piv.nlist(nnode) = 0;
          geom.piv.piv(nnode) = max(geom.piv.piv)+1;
        end % if( geom.piv.nlist(n(1))~=0 & geom.piv.nlist(n(2))~=0 )
        
      else % if( any(e==-1) )
        geom.piv.nlist(nnode) = 0;
        geom.piv.piv(nnode) = max(geom.piv.piv)+1;
      end %if( any(e==-1) )
    
    
    end % for l=1:geom.Nobj.N_edge

end