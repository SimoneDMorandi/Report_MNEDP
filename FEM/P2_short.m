
nnode = geom.Nobj.N_node;

for l=1:geom.Nobj.N_edge   %Processo lato per lato

  n(1)=geom.obj.E(l,1); % edge starting node
  n(2)=geom.obj.E(l,2); % edge ending node
  e(1)=geom.obj.E(l,3); % side element
  e(2)=geom.obj.E(l,4); % side element
  
  nnode = nnode + 1;

  geom.obj.P(nnode,:) = (geom.obj.P(n(1),:)+...
					geom.obj.P(n(2),:))/2; % P of the edge midpoint

  geom.obj.E(l,5)=nnode; % to connect the edge with its midpoint
  
  idx = [1 2 3];

  for el=e   % Ricerca del punto medio, fatto in modo ottimizzato
    
    if(el ~= -1)  % Se l elemento esiste
      acc = 0;
      acc = idx * ( geom.obj.T(el,1:3)==n(1) )'; %Guardo dov'è il vertice n-esimo nella matrice non estesa, resistuisce = [1,0,0] a seconda se vero o falso.
      acc = acc + idx * ( geom.obj.T(el,1:3)==n(2) )'; % Stessa cosa per il vertice di fine, poi sommo 

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

  VertexValue = [0 0];
%  Vertex = [0 0];
  D = [0 0];
  InputVertexValue=[0 0];
  
% Vedo se è un lato di D o N
% Se ai due estremi ho due nodi di D, devo distinguere due casi
	% Se in mezzo ho un nodo N, deduco che tutto il lato è di N e associo il gdl corrispondente.
	% Se il lato non è stato modificato, devo saperlo e dedurre dalle info iniziali il lato è N o D.

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


