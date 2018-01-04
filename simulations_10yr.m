%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina que carrega os dados do clima e clima_10p, faz as medias zonais  %
% e trimestrais e salva em um arquivo .mat                                %
% Autor: Livia Sancho                                                     %
% Data: Marco/2017                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clear all; close all; clc


% %LINUX
dir_trabalho='/media/SAMSUNG/Artigo_CESM_Tese_rafa/figuras';
mydados_clima='/media/SAMSUNG/cesm1_2_1/BCN_2/resultados_cam';
dir_dados = '/media/SAMSUNG/Artigo_CESM_Tese_rafa/dados';


VARIAVEL = {'PSL' 'TREFHT' 'PRECL' 'PRECC' 'LHFLX' 'CLDTOT' 'FLUT' 'OMEGA' 'U' 'V'};
variavel = {'psl' 'ts' 'precl' 'precc' 'lhflx' 'cldtot' 'flut' 'omega' 'u' 'v'}; %eh minusculo para dar o squeeze
RODADA = {'BCN_2' 'BCN_constsol_10p'};
COORDENADAS = {'lev' 'lon' 'lat'};
TRIMESTRE = {'JFM' 'AMJ' 'JAS' 'OND'};
NIVEL = {'13' '23'}; % referente a 200 e 850 hPa, respectivamente
nivel = {'200' '850'};


%% Clima

eval(['cd ',mydados_clima,'']) %entra no diretorio dos resultados do BCN_2

% Carrega cada variavel de cada arquivo nos 10 anos
for quantvariavel=1:length(VARIAVEL)
    for a=1:10 % anos
        for m=1:12 % meses
            ano=num2str(a, '%04i'); % identifica e transforma o numero em str
            mes=num2str(m, '%02i');
            %             Carrega os arquivos 2D
            if quantvariavel < 9
                eval([RODADA{1},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,' = ncread(''BCN_2.cam.h0.',ano,'-',mes,'.nc'',''',VARIAVEL{quantvariavel},''');'])
                %             Carrega os arquivos 3D
            else
                eval([RODADA{1},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,' = ncread(''BCN_2.cam.h0.',ano,'-',mes,'.nc'',''',VARIAVEL{quantvariavel},''');'])
                for n = 1:length(NIVEL);
                    eval([RODADA{1},'_',VARIAVEL{quantvariavel},'_',nivel{n},'_',ano,'_',mes,' = squeeze(',RODADA{1},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,'(:,:,',NIVEL{n},'));']) %minusculo: por causa do SQUEEZE
                end
            end
        end
    end
end


% Carrega as variaveis de {COORDENADAS}
for cont = 1:length(COORDENADAS)
    eval([COORDENADAS{cont},' = ncread(''BCN_2.cam.h0.0001-01.nc'',''',COORDENADAS{cont},''');'])
end


% Calcula as medias de acordo com a latitude
for cont=1:length(VARIAVEL);
    for a = 1:10
        for m = 1:12
            ano=num2str(a, '%04i');
            mes=num2str(m, '%02i');
            for i = 1:48 %pontos de lat
                if cont < 9;
                    eval(['x = squeeze(',RODADA{1},'_',VARIAVEL{cont},'_',ano,'_',mes,'(:,i));'])
                    eval([RODADA{1},'_',variavel{cont},'_',ano,'_',mes,'(i) = mean(x);'])
                else
                    for n = 1:length(NIVEL);
                        eval(['x = squeeze(',RODADA{1},'_',VARIAVEL{cont},'_',nivel{n},'_',ano,'_',mes,'(:,i));'])
                        eval([RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_',mes,'(i) = mean(x);'])
                    end
                end
            end
        end
    end
    eval(['clear ',RODADA{1},'_',VARIAVEL{cont},'_*;'])
end


% Calcula as medias trimestrais (JFM, AMJ, JAS, OND)
for a = 1:10
    for m = 1:3:12
        
        if m == 1
            p = 1;
        elseif m == 4
            p = 2;
        elseif m == 7
            p = 3;
        elseif m == 10
            p = 4;
        end
        
        ano=num2str(a, '%04i');
        mes=num2str(m, '%02i');
        
        for cont = 1:length(variavel) %variaveis com squeeze
            if cont < 9;
                if p < 4
                    eval(['iu = cat(1,',RODADA{1},'_',variavel{cont},'_',ano,'_',mes,',',RODADA{1},'_',variavel{cont},...
                        '_',ano,'_0',num2str(m)+1,',',RODADA{1},'_',variavel{cont},'_',ano,'_0',num2str(m)+2,');'])
                    eval([RODADA{1},'_',variavel{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                elseif p == 4
                    eval(['iu = cat(1,',RODADA{1},'_',variavel{cont},'_',ano,'_10,',RODADA{1},'_',variavel{cont},...
                        '_',ano,'_11,',RODADA{1},'_',variavel{cont},'_',ano,'_12);'])
                    eval([RODADA{1},'_',variavel{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                end
            else
                for n = 1:length(NIVEL);
                    if p < 4
                        eval(['iu = cat(1,',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_',mes,',',RODADA{1},'_',variavel{cont},...
                            '_',nivel{n},'_',ano,'_0',num2str(m)+1,',',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_0',num2str(m)+2,');'])
                        eval([RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                    elseif p == 4
                        eval(['iu = cat(1,',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_10,',RODADA{1},'_',variavel{cont},...
                            '_',nivel{n},'_',ano,'_11,',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_12);'])
                        eval([RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                    end
                end
            end
        end
    end
end


%calcula as medias trimestrais das rodadas
for cont = 1:length(variavel) %variaveis com squeeze
    for t = 1:length(TRIMESTRE);
        if cont < 9
            for a = 1:10;
                ano=num2str(a, '%04i');
                eval([RODADA{1},'_',variavel{cont},'_',TRIMESTRE{t},'(a,:) = ',RODADA{1},'_',variavel{cont},'_',ano,'_',TRIMESTRE{t},';'])
            end
            eval(['MEDIA','_',RODADA{1},'_',variavel{cont},'_',TRIMESTRE{t},' = nanmean(',RODADA{1},'_',variavel{cont},'_',TRIMESTRE{t},',1);'])
        else
            for n = 1:length(NIVEL);
                for a = 1:10;
                    ano=num2str(a, '%04i');
                    eval([RODADA{1},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},'(a,:) = ',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{t},';'])
                end
                eval(['MEDIA','_',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},' = nanmean(',RODADA{1},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},',1);'])
            end
        end
    end
end


% Limpa as variaveis que nao serao mais utilizadas
for v = 1:length(variavel);
    eval(['clear BCN_2_',variavel{v},'_0*;'])
end



% Salva os dados gerados em .mat

eval(['cd ',dir_dados,';'])

save clima_medias_trimestrais.mat lat lon lev

for v = 1:length(variavel);
    for t = 1:length(TRIMESTRE);
        if v < 9
            eval(['save (''clima_medias_trimestrais.mat'', ''MEDIA_BCN_2_',variavel{v},'_',TRIMESTRE{t},''',''-append'');'])
        else
            for n = 1:length(NIVEL);
                eval(['save (''clima_medias_trimestrais.mat'', ''MEDIA_BCN_2_',variavel{v},'_',nivel{n},'_',TRIMESTRE{t},''',''-append'');'])
            end
        end
    end
end



%% Clima_10p

clear all,close all,clc


% LINUX
dir_trabalho='/media/SAMSUNG/Artigo_CESM_Tese_rafa/figuras';
mydados_clima_10p='/media/SAMSUNG/cesm1_2_1/BCN_constsol_10p/run';
dir_dados = '/media/SAMSUNG/Artigo_CESM_Tese_rafa/dados';


VARIAVEL = {'PSL' 'TS' 'PRECL' 'PRECC' 'LHFLX' 'CLDTOT' 'FLUT' 'OMEGA' 'U' 'V'};
variavel = {'psl' 'ts' 'precl' 'precc' 'lhflx' 'cldtot' 'flut' 'omega' 'u' 'v'}; %eh minusculo para dar o squeeze
RODADA = {'BCN_2' 'BCN_constsol_10p'};
COORDENADAS = {'lev' 'lon' 'lat'};
TRIMESTRE = {'JFM' 'AMJ' 'JAS' 'OND'};
NIVEL = {'13' '23'}; % referente a 200 e 850 hPa, respectivamente
nivel = {'200' '850'};

eval(['cd ',mydados_clima_10p,'']) %entra no diretorio dos resultados do clima_10p

% Carrega cada variavel de cada arquivo nos 10 anos
for quantvariavel=1:length(VARIAVEL)
    for a=1:10 % anos
        for m=1:12 % meses
            ano=num2str(a, '%04i'); % identifica e transforma o numero em str
            mes=num2str(m, '%02i');
            %             Carrega os arquivos 2D
            if quantvariavel < 9
                eval([RODADA{2},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,' = ncread(''BCN_constsol_10p.cam.h0.',ano,'-',mes,'.nc'',''',VARIAVEL{quantvariavel},''');'])
                %             Carrega os arquivos 3D
            else
                eval([RODADA{2},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,' = ncread(''BCN_constsol_10p.cam.h0.',ano,'-',mes,'.nc'',''',VARIAVEL{quantvariavel},''');'])
                for n = 1:length(NIVEL);
                    eval([RODADA{2},'_',VARIAVEL{quantvariavel},'_',nivel{n},'_',ano,'_',mes,' = squeeze(',RODADA{2},'_',VARIAVEL{quantvariavel},'_',ano,'_',mes,'(:,:,',NIVEL{n},'));']) %minusculo: por causa do SQUEEZE
                end
            end
        end
    end
end


% Carrega as variaveis de {COORDENADAS}
for cont = 1:length(COORDENADAS)
    eval([COORDENADAS{cont},' = ncread(''BCN_constsol_10p.cam.h0.0001-01.nc'',''',COORDENADAS{cont},''');'])
end


% Calcula as medias de acordo com a latitude
for cont=1:length(VARIAVEL)
    for a = 1:10
        for m = 1:12
            ano=num2str(a, '%04i');
            mes=num2str(m, '%02i');
            for i = 1:48 %pontos de lat
                if cont < 9;
                    eval(['x = squeeze(',RODADA{2},'_',VARIAVEL{cont},'_',ano,'_',mes,'(:,i));'])
                    eval([RODADA{2},'_',variavel{cont},'_',ano,'_',mes,'(i) = mean(x);'])
                else
                    for n = 1:length(NIVEL);
                        eval(['x = squeeze(',RODADA{2},'_',VARIAVEL{cont},'_',nivel{n},'_',ano,'_',mes,'(:,i));'])
                        eval([RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_',mes,'(i) = mean(x);'])
                    end
                end
            end
        end
    end
end


% Calcula as medias trimestrais (JFM, AMJ, JAS, OND)
for a = 1:10
    for m = 1:3:12
        
        if m == 1
            p = 1;
        elseif m == 4
            p = 2;
        elseif m == 7
            p = 3;
        elseif m == 10
            p = 4;
        end
        
        ano=num2str(a, '%04i');
        mes=num2str(m, '%02i');
        
        for cont = 1:length(variavel) %variaveis com squeeze
            if cont < 9;
                if p < 4
                    eval(['iu = cat(1,',RODADA{2},'_',variavel{cont},'_',ano,'_',mes,',',RODADA{2},'_',variavel{cont},...
                        '_',ano,'_0',num2str(m)+1,',',RODADA{2},'_',variavel{cont},'_',ano,'_0',num2str(m)+2,');'])
                    eval([RODADA{2},'_',variavel{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                elseif p == 4
                    eval(['iu = cat(1,',RODADA{2},'_',variavel{cont},'_',ano,'_10,',RODADA{2},'_',variavel{cont},...
                        '_',ano,'_11,',RODADA{2},'_',variavel{cont},'_',ano,'_12);'])
                    eval([RODADA{2},'_',variavel{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                end
            else
                for n = 1:length(NIVEL);
                    if p < 4
                        eval(['iu = cat(1,',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_',mes,',',RODADA{2},'_',variavel{cont},...
                            '_',nivel{n},'_',ano,'_0',num2str(m)+1,',',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_0',num2str(m)+2,');'])
                        eval([RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                    elseif p == 4
                        eval(['iu = cat(1,',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_10,',RODADA{2},'_',variavel{cont},...
                            '_',nivel{n},'_',ano,'_11,',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_12);'])
                        eval([RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                    end
                end
            end
        end
    end
end


%calcula as medias trimestrais das rodadas
for cont = 1:length(variavel) %variaveis com squeeze
    for t = 1:length(TRIMESTRE);
        
        if cont < 9
            for a = 1:10;
                ano=num2str(a, '%04i');
                eval([RODADA{2},'_',variavel{cont},'_',TRIMESTRE{t},'(a,:) = ',RODADA{2},'_',variavel{cont},'_',ano,'_',TRIMESTRE{t},';'])
            end
            eval(['MEDIA','_',RODADA{2},'_',variavel{cont},'_',TRIMESTRE{t},' = nanmean(',RODADA{2},'_',variavel{cont},'_',TRIMESTRE{t},',1);'])
        else
            for n = 1:length(NIVEL);
                for a = 1:10;
                    ano=num2str(a, '%04i');
                    eval([RODADA{2},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},'(a,:) = ',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',ano,'_',TRIMESTRE{t},';'])
                end
                eval(['MEDIA','_',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},' = nanmean(',RODADA{2},'_',variavel{cont},'_',nivel{n},'_',TRIMESTRE{t},',1);'])
            end
        end
    end
end


% Limpa as variaveis que nao serao mais utilizadas
for v = 1:length(variavel);
    eval(['clear BCN_constsol_10p_',variavel{v},'_0*;'])
end


% Salva os dados gerados em .mat

eval(['cd ',dir_dados,';'])

save clima_10p_medias_trimestrais.mat lat

for v = 1:length(variavel);
    for t = 1:length(TRIMESTRE);
        if v < 9
            eval(['save (''clima_10p_medias_trimestrais.mat'', ''MEDIA_BCN_constsol_10p_',variavel{v},'_',TRIMESTRE{t},''',''-append'');'])
        else
            for n = 1:length(NIVEL);
                eval(['save (''clima_10p_medias_trimestrais.mat'', ''MEDIA_BCN_constsol_10p_',variavel{v},'_',nivel{n},'_',TRIMESTRE{t},''',''-append'');'])
            end
        end
    end
end


disp('Chegou ao final de separa_dados_clima_clima_10p.m com sucesso!')