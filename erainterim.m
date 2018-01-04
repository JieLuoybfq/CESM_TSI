%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina que carrega os dados do ERA-INTERIM, faz as medias zonais e      %
% trimestrais e salva em um arquivo .mat                                  %
% Autor: Livia Sancho                                                     %
% Data: Marco/2017                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clear all; close all; clc

%WINDOWS
dir_trabalho='D:\Artigo_CESM_Tese_rafa\figuras';
dir_dados = 'D:\Artigo_CESM_Tese_rafa\dados';
dir_erainterim='D:\Artigo_CESM_Tese_rafa\era_interim\netcdf';


% Carrega os dados de slp do ERA
eval(['cd ',dir_erainterim,'']) %entra no diretorio da Reanalise
% dado de lat
ERAINTERIM_LAT = ncread('ei.moda.an.sfc.regn128sc.2000010100.nc','g4_lat_0');

% dado de lat
ERAINTERIM_LON = ncread('ei.moda.an.sfc.regn128sc.2000010100.nc','g4_lon_1');

% (Mean sea level pressure, 2 metre temperature)
VARIAVEL_ERAINTERIM = {'MSL_GDS4_SFC_123' '2T_GDS4_SFC_123'};
variavel_erainterim = {'msl_gds4_sfc_123' '2t_gds4_sfc_123'}; %eh minusculo para dar o squeeze
RODADA_ERAINTERIM = {'ERAINTERIM'};
TRIMESTRE = {'JFM' 'AMJ' 'JAS' 'OND'};

% Carrega os arquivos do ERA e faz a media na longitude
for quantvariavel_era=1:length(VARIAVEL_ERAINTERIM)
    for a=2000:2009 % anos
        for m=1:12 % meses
            ano=num2str(a, '%04i');
            mes=num2str(m, '%02i');
            eval([RODADA_ERAINTERIM{1},'_',VARIAVEL_ERAINTERIM{quantvariavel_era},'_',ano,'_',mes,' = ncread(''ei.moda.an.sfc.regn128sc.',ano,'',mes,'0100.nc'',''',VARIAVEL_ERAINTERIM{quantvariavel_era},''');'])
            for i = 1:256 %pontos de lat
                eval(['x = squeeze(',RODADA_ERAINTERIM{1},'_',VARIAVEL_ERAINTERIM{quantvariavel_era},'_',ano,'_',mes,'(:,i));'])
                eval([RODADA_ERAINTERIM{1},'_',variavel_erainterim{quantvariavel_era},'_',ano,'_',mes,'(i) = nanmean(x);'])
            end
        end
    end
end

clear ERAINTERIM_*GDS4*


% Calcula as medias trimestrais (JFM, AMJ, JAS, OND)
for a = 2000:2009 % anos
    for m = 1:3:12 % meses de 3 em 3
        
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
        for cont = 1:length(variavel_erainterim) %variaveis com squeeze
            if p < 4
                eval(['iu = cat(1,',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_',mes,',',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},...
                    '_',ano,'_0',num2str(m)+1,',',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_0',num2str(m)+2,');'])
                eval([RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                clear iu
            elseif p == 4
                eval(['iu = cat(1,',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_10,',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},...
                    '_',ano,'_11,',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_12);'])
                eval([RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_',TRIMESTRE{p},' = mean(iu,1);'])
                clear iu
            end
        end
    end
end



%calcula as medias trimestrais das rodadas
for cont = 1:length(variavel_erainterim) %variaveis com squeeze
    for t = 1:length(TRIMESTRE);
        for a = 1:10;
            ano=num2str(a-1+2000, '%04i');
            eval([RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',TRIMESTRE{t},'(a,:) = ',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',ano,'_',TRIMESTRE{t},';'])
        end
        eval(['MEDIA','_',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',TRIMESTRE{t},' = nanmean(',RODADA_ERAINTERIM{1},'_',variavel_erainterim{cont},'_',TRIMESTRE{t},',1);'])
    end
end



% Salva os dados gerados em .mat
eval(['cd ',dir_dados,';'])


TRIMESTRE = {'JFM' 'AMJ' 'JAS' 'OND'};

save era_medias_trimestrais.mat ERAINTERIM_LAT

for v = 1:length(variavel_erainterim);
    for t = 1:length(TRIMESTRE);
    eval(['save (''era_medias_trimestrais.mat'', ''MEDIA_ERAINTERIM_',variavel_erainterim{v},'_',TRIMESTRE{t},''',''-append'');'])
    end
end

disp('Chegou ao final de separa_dados_era.m com sucesso!')


