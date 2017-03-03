function [controls,markers,control_field,marker_field]=populate_search(method)
%----
% generate key for strain normalization
%----
    switch method
        case 'RNAP_marker'
            controls.rpoBC={'rpoBC{cat}'};
            controls.dksA={'rpoBC{cat};dksA'};
            controls.rpoD ={'rpoD{cat}'};
            controls.fis  ={'fis{aph}'};
            controls.nusAaph ={'nusA{aph}'};
            controls.nusAcat ={'nusA{cat}'};
            controls.nusGaph ={'nusG{aph}'};
            controls.nusGcat ={'nusG{cat}'};
            controls.rhoaph  ={'rho{aph}'};
            controls.rhocat  ={'rho{cat}'};
            controls.nusE ={'nusE{cat}'};
            controls.rpoC ={'rpoC{cat}'};
                      
            markers.rpoBC ={'ALS2*3'};
            markers.dksA ={'ALS2*3; aph(KEIO:dksA)'};
            markers.rpoD ={'ALS4*5','ALS13*14','ALS49*50','aph(AH:rpoD)'};
            markers.fis  ={'ALS90*91','ALS92*93'};
            markers.nusAaph ={'aph(JP:nusA)'};
            markers.nusAcat ={'cat(JP:nusA)'};
            markers.nusGaph ={'aph(JP:nusG)'};
            markers.nusGcat ={'cat(JP:nusG)'};
            markers.rhoaph  ={'aph(JP:rho)'};
            markers.rhocat  ={'cat(JP:rho)'};
            markers.nusE ={'cat(JP:nusE)'};
            markers.rpoC ={'cat(WR:rpoC)'};
            
            control_field='mut';
            marker_field='mrk';
            
        case 'marker'
            controls.rpoD ={'rpoD{cat}'};
            controls.fis  ={'fis{aph}'};
            controls.nusAaph ={'nusA{aph}'};
            controls.nusAcat ={'nusA{cat}'};
            controls.nusGaph ={'nusG{aph}'};
            controls.nusGcat ={'nusG{cat}'};
            controls.rhoaph  ={'rho{aph}'};
            controls.rhocat  ={'rho{cat}'};
            controls.nusE ={'nusE{cat}'};
            controls.rpoC ={'rpoC{cat}'};
                      
            markers.rpoD ={'ALS4*5','ALS13*14','ALS49*50','aph(AH:rpoD)'};
            markers.fis  ={'ALS90*91','ALS92*93'};
            markers.nusAaph ={'aph(JP:nusA)'};
            markers.nusAcat ={'cat(JP:nusA)'};
            markers.nusGaph ={'aph(JP:nusG)'};
            markers.nusGcat ={'cat(JP:nusG)'};
            markers.rhoaph  ={'aph(JP:rho)'};
            markers.rhocat  ={'cat(JP:rho)'};
            markers.nusE ={'cat(JP:nusE)'};
            markers.rpoC ={'cat(WR:rpoC)'};
            
            control_field='mut';
            marker_field='mrk';
            
        case 'RNAPgroup'
            
            controls.rpoBC={'ALS2*3'};
            markers.rpoBC={'ALS2*3','ALS2*3; aph(KEIO:dksA)'};
            
            control_field='mrk';
            marker_field='mrk';
            
        case 'RNAPcontrol'
            controls.rpoBC={'rpoBC{cat}'};
            
             markers.rpoBC={'rpoBC{cat}','rpoB{D1297E}','rpoB{I572T,S574F}',...
                            'rpoB{D516N}','rpoB{S574F}','rpoB{S574Y}',...
                            'rpoC{R88C}','rpoC{N690C}','rpoC{G729D}',...
                            'rpoC{V1141S}','rpoB{I1112S}','rpoC{D643G}',...
                            'rpoC{K650T}','rpoB{225-GG-343}','rpoB{D516V}',...
                            'rpoB{N518D}','rpoB{Q513L}','rpoB{del(507-511),ins(V),P552L}',...
                            'rpoB{S531F}','rpoB{V146F}','rpoB{Q513P}',...
                            'rpoB{V146G}','rpoB{I572F}'};
           
           control_field='mut';
           marker_field='mut';
           
        case 'pAsensitive'
            controls.rpoBC={'rpoB{G181V}','rpoC{R1330S}','rpoC{D264Y,319ins(TKRPL)}',...
                            'rpoC{D264Y}','rpoB{L533H}','rpoB{Q148R}',...
                            'rpoB{S512(F|Y),I572Q}','rpoB{P564L}','rpoC{L796T}',...
                            'rpoC{N458A}','rpoC{P1022L}','rpoB{P1081L}',...
                            'rpoC{K445P}','rpoB{937-GGG-1040}','rpoB{D446H}'};
            
             markers.rpoBC={'rpoB{G181V}','rpoC{R1330S}','rpoC{D264Y,319ins(TKRPL)}',...
                            'rpoC{D264Y}','rpoB{L533H}','rpoB{Q148R}',...
                            'rpoB{S512(F|Y),I572Q}','rpoB{P564L}','rpoC{L796T}',...
                            'rpoC{N458A}','rpoC{P1022L}','rpoB{P1081L}',...
                            'rpoC{K445P}','rpoB{937-GGG-1040}','rpoB{D446H}'};           
                         
            control_field='mut';
            marker_field='mut';
        
        case 'pAresistant'
            controls.rpoBC={'rpoB{Q148P}','rpoB{R454H}','rpoC{G1354C}',...
                            'rpoC{G333C}','rpoC{L1314Q}','rpoC{del(215-220)}',...
                            'rpoB{H447P}','rpoB{H551P}','rpoB{G570C}',...
                            'rpoB{T563P}','rpoB{del(532)}','rpoB{L533P}',...
                            'rpoB{L571Q}','rpoB{P153L}','rpoB{Y395D}',...
                            'rpoC{del(1338-1341)}','rpoC{S1324L}','rpoB{R451C}',...
                            'rpoB{del(442),D443S,D444M,D446S,H447P,L448S,R451C}',...
                            'rpoC{D348Y}','rpoB{G1260D}'};
             
             markers.rpoBC={'rpoB{Q148P}','rpoB{R454H}','rpoC{G1354C}',...
                            'rpoC{G333C}','rpoC{L1314Q}','rpoC{del(215-220)}',...
                            'rpoB{H447P}','rpoB{H551P}','rpoB{G570C}',...
                            'rpoB{T563P}','rpoB{del(532)}','rpoB{L533P}',...
                            'rpoB{L571Q}','rpoB{P153L}','rpoB{Y395D}',...
                            'rpoC{del(1338-1341)}','rpoC{S1324L}','rpoB{R451C}',...
                            'rpoB{del(442),D443S,D444M,D446S,H447P,L448S,R451C}',...
                            'rpoC{D348Y}','rpoB{G1260D}'};
            control_field='mut';
            marker_field='mut';
    end
    
  
end