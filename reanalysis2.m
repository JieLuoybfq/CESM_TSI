%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina que carrega os dados da Reanalise, faz as medias zonais e        %
% trimestrais e salva em um arquivo .mat                                  %
% Autor: Livia Sancho                                                     %
% Data: Marco/2017                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

clear all; close all; clc

%WINDOWS
dir_trabalho='D:\Artigo_CESM_Tese_rafa\figuras';
dir_dados = 'D:\Artigo_CESM_Tese_rafa\dados';
dir_reanalise='D:\Artigo_CESM_Tese_rafa\reanalise';

% Carrega os dados da REANALISE
eval(['cd ',dir_reanalise,'']) %entra no diretorio da Reanalise

nome_arq = {'slp' 'temp'}; % variavel do nome do arquivo
var = {'mslp' 'air'}; % variavel de dentro do arquivo
VAR = {'SLP' 'TEMP'}; % nome da variavel que sera criada

% Carrega os arquivos da REANALISE
for v = 1:length(var);
    eval(['REANALISE_',VAR{v},' = ncread(''',nome_arq{v},'_reanalisys2_jan2000-dec2010.nc'',''',var{v},''');'])
    if v == 2
        eval(['REANALISE_',VAR{v},' = squeeze(REANALISE_',VAR{v},');'])
    end
%     Tira os ultimos 12 meses do tempo pois o arquivo de leitura e de
%     2000-2010 e nao de 2000-2009
    eval(['REANALISE_',VAR{v},' = REANALISE_',VAR{v},'(:,:,1:end-12);'])
end



% dado de lat
REANALISE_LAT = ncread('slp_reanalisys2_jan2000-dec2009.nc','lat');
% dado de lon
REANALISE_LON = ncread('slp_reanalisys2_jan2000-dec2009.nc','lon');
% dado de tempo
REANALISE_TIME = ncread('slp_reanalisys2_jan2000-dec2009.nc','time');
REANALISE_TIME = REANALISE_TIME/24+datenum(1800,1,1); % Reanalise fornece o tempo (horas) iniciando em 01/01/1800.
REANALISE_TEMPO = datestr(REANALISE_TIME);


% Calcula as medias na longitude
% slp_media_lat = squeeze(nanmean(REANALISE_SLP,1));
for v = 1:length(var);
    eval([nome_arq{v},'_media_lat = squeeze(nanmean(REANALISE_',VAR{v},',1));'])
end

% Adquire meses do vetor de tempo
reanalise_mes = str2num(datestr(REANALISE_TIME,'mm'));


for v = 1:length(var);
%     Janeiro, Fevereiro e Marco (JFM)
    
%     identifica os meses de janeiro ou fevereiro ou marco e joga para a
%     variavel idx_mes
    idx_mes = find(reanalise_mes==1 | reanalise_mes==2 | reanalise_mes==3);
    
%     gera uma variavel *_jfm1 com todas as latitudes e os meses
%     identificados no comando anterior
    eval([nome_arq{v},'_jfm1 = ',nome_arq{v},'_media_lat(:,idx_mes);'])
    clear idx
%     Faz a media dos meses selecionados
    eval(['MEDIA_REANALISE_',nome_arq{v},'_JFM = nanmean(',nome_arq{v},'_jfm1,2);'])
    
%     Abril, Maio e Junho (AMJ)
    idx_mes = find(reanalise_mes==4 | reanalise_mes==5 | reanalise_mes==6);
    eval([nome_arq{v},'_amj1 = ',nome_arq{v},'_media_lat(:,idx_mes);'])
    clear idx
    MEDIA_REANALISE_slp_AMJ = nanmean(slp_amj1,2);
    eval(['MEDIA_REANALISE_',nome_arq{v},'_AMJ = nanmean(',nome_arq{v},'_amj1,2);'])
    
%     Julho, Agosto e Setembro (JAS)
    idx_mes = find(reanalise_mes==7| reanalise_mes==8 | reanalise_mes==9);
    eval([nome_arq{v},'_jas1 = ',nome_arq{v},'_media_lat(:,idx_mes);'])
    clear idx
    eval(['MEDIA_REANALISE_',nome_arq{v},'_JAS = nanmean(',nome_arq{v},'_jas1,2);'])
    
%     Outubro, Novembro e Dezembro (OND)
    idx_mes = find(reanalise_mes==10| reanalise_mes==11 | reanalise_mes==12);
    eval([nome_arq{v},'_ond1 = ',nome_arq{v},'_media_lat(:,idx_mes);'])
    clear idx
    MEDIA_REANALISE_slp_OND = nanmean(slp_ond1,2);
    eval(['MEDIA_REANALISE_',nome_arq{v},'_OND = nanmean(',nome_arq{v},'_ond1,2);'])
    
end


% Salva os dados gerados em .mat
eval(['cd ',dir_dados,';'])

TRIMESTRE = {'JFM' 'AMJ' 'JAS' 'OND'};

save reanalise_medias_trimestrais.mat REANALISE_LAT

for v = 1:length(var);
    for t = 1:length(TRIMESTRE);
        eval(['save (''reanalise_medias_trimestrais.mat'', ''MEDIA_REANALISE_',nome_arq{v},'_',TRIMESTRE{t},''',''-append'');'])
    end
end


disp('Chegou ao final de separa_dados_reanalise.m com sucesso!')
