%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotina que plota as figuras que serao usadas no artigo                  % 
% Autor: Livia Sancho                                                     %
% Data: Maio/2017                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

start

clear all; close all; clc



% LINUX
dir_trabalho='/media/SAMSUNG/Artigo_CESM_Tese_rafa/figuras';
dir_dados = '/media/SAMSUNG/Artigo_CESM_Tese_rafa/dados';
addpath(genpath('/home/numa2/Documents/MATLAB/toolbox/shadedErrorBar'));
addpath(genpath('/home/numa2/Documents/MATLAB/toolbox/mtit'));


base_dados = {'clima' 'clima_10p' 'era' 'reanalise'};


% Carrega os dados
eval(['cd ',dir_dados,''])

for b = 1:length(base_dados);
    eval(['load ',base_dados{b},'_medias_trimestrais.mat'])
%     clear *_OND* *_AMJ*
end


clear *_cldtot_* *_flut_* *_lhflx_* *_precc_* *_precl_* *_u_* *_v_* *_omega_*

era = {'msl' '2t'};
reanalise = {'slp' 'temp'};
cesm = {'psl' 'ts'};
trimestre = {'JFM' 'AMJ' 'JAS' 'OND'};
ylbl1 = {'Mean Sea' 'Surface Air'};
ylbl2 = {'Level Pressure (hPa)' 'Temperature (°C)'};
% ttl = {'Mean Sea Level Pressure' 'Surface Air Temperature'};


letras = {'(a)' '' '(b)' ''};

%% Figuras clima, reanalise e era em um mesmo subplot

letras_sub = {'(a)' '(b)' '(c)' '(d)'};
sub = {'1' '2' '3' '4'};


figure
for v = 1:length(cesm);
    for t = 1:2:length(trimestre);
        
        if v == 1 && t == 1;
            h1 = subplot(2,2,1)
        elseif v == 1 && t == 3;
            h2 = subplot(2,2,2)
        elseif v == 2 && t == 1;
            h3 = subplot(2,2,3)
        elseif v == 2 && t == 3;
            h4 = subplot(2,2,4)
        end
        
        if v == 2
            eval(['p1=plot(lat,MEDIA_BCN_2_',cesm{v},'_',trimestre{t},'-273.15,''b'',''linewidth'',1.5);'])
            hold on
            eval(['p3=plot(REANALISE_LAT,MEDIA_REANALISE_',reanalise{v},'_',trimestre{t},'-273.15,''k--'',''linewidth'',1.5);'])
            eval(['p4=plot(ERAINTERIM_LAT,MEDIA_ERAINTERIM_',era{v},'_gds4_sfc_123_',trimestre{t},'-273.15,''g--'',''linewidth'',1.5);'])
            
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            axis([-88 88 -60 40])
            
        else
            eval(['p1=plot(lat,MEDIA_BCN_2_',cesm{v},'_',trimestre{t},'/100,''b'',''linewidth'',1.5);'])
            hold on
            eval(['p3=plot(REANALISE_LAT,MEDIA_REANALISE_',reanalise{v},'_',trimestre{t},'/100,''k--'',''linewidth'',1.5);'])
            eval(['p4=plot(ERAINTERIM_LAT,MEDIA_ERAINTERIM_',era{v},'_gds4_sfc_123_',trimestre{t},'/100,''g--'',''linewidth'',1.5);'])
            
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            axis([-88 88 980 1025])
            
        end
        
        plot([-88 88],[0 0],'k--')
        
        %         eval(['text(0.05,0.9,''',letras_sub{s},''',''Units'',''normalized'');'])
        
        grid on
        
        if t == 1
            eval(['ylabel({''',ylbl1{v},''';''',ylbl2{v},'''},''fontsize'',14)'])
        end
        
        if v == 2 && t == 3
%             hh = mtit('bla');
%             set(hh, 'visible', 'off')
%             xlh = xlabel(hh.ah, 'Latitude');
%             set(xlh, 'Visible', 'On','fontsize',14)
            xlabel('Latitude','fontsize',14)
        end
        
      
        
                if t == 3 && v == 2;
                    leg = legend([p1,p3,p4],{'<CONTROL>','Reanalisys2','Era-Interim'},'location','bestoutside');%,'location','southoutside','orientation','horizontal');
                end
        
        pbaspect([1 .3 1])
    end
end


tightfig

cd ../figuras

eval(['namest = ''',cesm{v},'_clima_reanalise_era_cesm_subplot_LS'';'])
saveas(gcf,namest,'fig')
print(namest,'-dpng','-r300');
close   


%% Figuras com desvio padrao e bias no mesmo grafico

clear all

% WINDOWS
% dir_trabalho = 'D:\Artigo_CESM_Tese_rafa\figuras';
% dir_dados = 'D:\Artigo_CESM_Tese_rafa\dados';

% LINUX
dir_trabalho='/media/SAMSUNG/Artigo_CESM_Tese_rafa/figuras';
dir_dados = '/media/SAMSUNG/Artigo_CESM_Tese_rafa/dados';


base_dados = {'clima' 'clima_10p'};
rodada = {'BCN_2' 'BCN_constsol_10p'};

cesm = {'psl' 'ts' 'lhflx' 'cldtot' 'precc' 'precl' 'u_200' 'u_850' 'v_200' 'v_850'};

ylbl = {'Mean Sea Level Pressure (hPa)' 'Surface Temperature (°C)' 'Surface Latent Heat Flux (Wm^-^2)' ...
    'Total Cloud Cover (%)' 'Precipitation (mm month^-^1)' 'Zonal Wind 200 hPa (m^-^s)' 'Zonal Wind 850 hPa (m^-^s)' ...
    'Meridional Wind 200 hPa (m^-^s)' 'Meridional Wind 850 hPa (m^-^s)'};
% ylbl2 = {'Pressure (hPa)' 'Temperature (°C)' 'Heat Flux (Wm^-^2)' 'Cloud Cover (%)' ...
%     'Wind (m^-^s)' 'Wind (m^-^s)' 'Precipitation (mm month^-^1)'};
ttl = {'Mean Sea Level Pressure' 'Surface Temperature' 'Surface Latent Heat Flux' ...
    'Total Cloud Cover' 'Precipitation' 'Zonal Wind 200 hPa' 'Zonal Wind 850 hPa' 'Meridional Wind 200 hPa' ...
    'Meridional Wind 850 hPa'};

trimestre = {'JFM' 'AMJ' 'JAS' 'OND'};


% Carrega os dados
eval(['cd ',dir_dados,''])

for b = 1:length(base_dados);
    eval(['load ',base_dados{b},'_medias_trimestrais.mat'])
%     clear *_OND* *_AMJ*
end

clear b *_omega_*


% Calcula a precipitacao total por mes por rodada por trimestre
for r = 1:length(rodada);
    for t = 1:length(trimestre)
        eval(['MEDIA_',rodada{r},'_prec_',trimestre{t},' = (MEDIA_',rodada{r},'_precc_',trimestre{t}, ...
            ' + MEDIA_',rodada{r},'_precl_',trimestre{t},')*3600*24*30*1000;'])
    end
end
clear r t





clear *REANALISE_* *ERAINTERIM_*

% Carrega os demais dados
for b = 1:length(base_dados); %clima e clima_10p
    eval(['load ',base_dados{b},'_medias_trimestrais_anual.mat'])
%     clear *_OND* *_AMJ*
end


cesm2 = {'psl' 'ts' 'lhflx' 'cldtot' 'prec' 'u_200' 'u_850' 'v_200' 'v_850'};


for r = 1:length(rodada); % rodadas - clima e clima_10p
    for v = 1:length(cesm2); % variaveis do cesm
        for t = 1:length(trimestre); % trimestres
            
            if v == 5;
                % Calcula a precipitacao total por mes por rodada por trimestre
                eval(['MEDIA_',rodada{r},'_prec_',trimestre{t},' = (MEDIA_',rodada{r},'_precc_',trimestre{t}, ...
                    ' + MEDIA_',rodada{r},'_precl_',trimestre{t},')*3600*24*30*1000;'])
                
                eval([rodada{r},'_prec_',trimestre{t},' = (',rodada{r},'_precc_',trimestre{t}, ...
                    ' + ',rodada{r},'_precl_',trimestre{t},')*3600*24*30*1000;'])
                
%                 cesm2(v) = {'prec'};
            end

            eval(['std_',rodada{r},'_',cesm2{v},'_',trimestre{t},' = std(',rodada{r},'_',cesm2{v},'_',trimestre{t},');'])
            
            eval([rodada{r},'_',cesm2{v},'_',trimestre{t},'1 = MEDIA_',rodada{r},'_',cesm2{v},'_',trimestre{t}, ...
                ' - std_',rodada{r},'_',cesm2{v},'_',trimestre{t},';'])
            
            eval([rodada{r},'_',cesm2{v},'_',trimestre{t},'2 = MEDIA_',rodada{r},'_',cesm2{v},'_',trimestre{t}, ...
                ' + std_',rodada{r},'_',cesm2{v},'_',trimestre{t},';'])
            
            eval(['MEDIA_',rodada{r},'_',cesm2{v},'_',trimestre{t},'_std = [',rodada{r},'_',cesm2{v},'_',trimestre{t}, ...
                '1;MEDIA_',rodada{r},'_',cesm2{v},'_',trimestre{t},';',rodada{r},'_',cesm2{v},'_',trimestre{t},'2];'])
        end
    end
end


axish1min = {'975' '-50' '0' '30' '0' '-15' '-10' '-4.5' '-2'};
axish1max = {'1030' '50' '165' '100' '300' '60' '20' '4.5' '2'};


axish2min = {'-10' '0' '-10' '-10' '-20' '-10' '-10' '-1' '-0.3'};
axish2max = {'10' '20' '20' '20' '40' '15' '10' '1' '0.3'};


 


letras = {'(a)' '' '(b)' ''};

left= 0.15;
bottom1=0.5;
bottom2=0.25;
width=0.8;
height1=0.45;% which is also bottom1-bottom2
height2=0.248;

ylbl2 = {'Mean Sea' 'Surface Air' 'Surface Latent' 'Total Cloud' 'Precipitation' 'Zonal Wind' 'Zonal Wind' ...
    'Meridional Wind' 'Meridional Wind'};

ylbl3 = {'Level Pressure (hPa)' 'Temperature (°C)' 'Heat Flux (Wm^-^2)' 'Cover (%)' '(mm month^-^1)' '200 hPa (m^-^s)' '850 hPa (m^-^s)' ...
    '200 hPa (m^-^s)' '850 hPa (m^-^s)'};

% stop

for v = 1:length(cesm2); % variaveis do cesm
    for t = 1:length(trimestre); % trismestres
        
        figure
        if v == 1;
            axes('Position',[left bottom1 width height1]);
            eval(['H(1) = shadedErrorBar(lat,MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'_std/100, {@mean, @(x) 1*std(x)}, {''-b'', ''LineWidth'', 1.5}, 0);'])
            hold on
            eval(['H(2) = shadedErrorBar(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'_std/100, {@mean, @(x) 1*std(x)}, {''-r'', ''LineWidth'', 1.5}, 0);'])
            plot([-88 88],[0 0],'k--')
            %             legend([H(1).mainLine, H(2).mainLine, H.patch, H.patch],'<CONTROL>','<TSI10p>','STD','STD');
            if t == 3;
                legend([H(1).mainLine, H(2).mainLine],'<CONTROL>','<TSI10p>','location','best');
            end
            %             pbaspect([1 .3 1])
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            eval(['axis([-88 88 ',axish1min{v},' ',axish1max{v},'])'])
%             eval(['title(''',ttl{v},' - ',trimestre{t},''',''fontsize'',11)'])
            eval(['h_label=ylabel({''',ylbl2{v},''';''',ylbl3{v},'''},''fontsize'',14,''visible'',''on'');'])
            set(gca, 'XTickLabel', [],'XTick',[])
            eval(['set(gca,''ytick'',[',axish1min{v},'+10:10:',axish1max{v},']);'])
            
            eval(['text(0.05,0.9,''',letras{t},''',''Units'',''normalized'');'])
            
            axes('Position',[left bottom2 width height2])
            eval(['plot(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'/100 - MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'/100,''k'',''linewidth'',1.5);'])
            hold on
            plot([-88 88],[0 0],'k--')
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            set(gca,'YAxisLocation','right');
            eval(['axis([-88 88 ',axish2min{v},' ',axish2max{v},'])'])
            ylabel('Difference','fontsize',14)
            xlabel('Latitude','fontsize',14)
            
            tightfig
            
        elseif v == 2;
            axes('Position',[left bottom1 width height1]);
            eval(['H(1) = shadedErrorBar(lat,MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'_std-273.15, {@mean, @(x) 1*std(x)}, {''-b'', ''LineWidth'', 1.5}, 0);'])
            hold on
            eval(['H(2) = shadedErrorBar(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'_std-273.15, {@mean, @(x) 1*std(x)}, {''-r'', ''LineWidth'', 1.5}, 0);'])
            plot([-88 88],[0 0],'k--')
            if t == 3;
                legend([H(1).mainLine, H(2).mainLine],'<CONTROL>','<TSI10p>','location','best');
            end
            %             pbaspect([1 .3 1])
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            eval(['axis([-88 88 ',axish1min{v},' ',axish1max{v},'])'])
%             eval(['title(''',ttl{v},' - ',trimestre{t},''',''fontsize'',11)'])
            eval(['h_label=ylabel({''',ylbl2{v},''';''',ylbl3{v},'''},''fontsize'',14,''visible'',''on'');'])
            set(gca, 'XTickLabel', [],'XTick',[])
            eval(['set(gca,''ytick'',[',axish1min{v},'+10:10:',axish1max{v},']);'])
            
            eval(['text(0.05,0.9,''',letras{t},''',''Units'',''normalized'');'])
            
            axes('Position',[left bottom2 width height2])
            eval(['plot(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'-273.15 - (MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'-273.15),''k'',''linewidth'',1.5);'])
            hold on
            plot([-88 88],[0 0],'k--')
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            set(gca,'YAxisLocation','right');
            eval(['axis([-88 88 ',axish2min{v},' ',axish2max{v},'])'])
            ylabel('Difference','fontsize',14)
            xlabel('Latitude','fontsize',14)
            
            tightfig
            
        elseif v == 4;
            axes('Position',[left bottom1 width height1]);
            eval(['H(1) = shadedErrorBar(lat,MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'_std*100, {@mean, @(x) 1*std(x)}, {''-b'', ''LineWidth'', 1.5}, 0);'])
            hold on
            eval(['H(2) = shadedErrorBar(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'_std*100, {@mean, @(x) 1*std(x)}, {''-r'', ''LineWidth'', 1.5}, 0);'])
            plot([-88 88],[0 0],'k--')
            if t == 3;
                legend([H(1).mainLine, H(2).mainLine],'<CONTROL>','<TSI10p>','location','best');
            end
            %             pbaspect([1 .3 1])
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            eval(['axis([-88 88 ',axish1min{v},' ',axish1max{v},'])'])
%             eval(['title(''',ttl{v},' - ',trimestre{t},''',''fontsize'',11)'])
            eval(['h_label=ylabel({''',ylbl2{v},''';''',ylbl3{v},'''},''fontsize'',14,''visible'',''on'');'])
            set(gca, 'XTickLabel', [],'XTick',[])
            eval(['set(gca,''ytick'',[',axish1min{v},'+10:10:',axish1max{v},']);'])
            
            eval(['text(0.05,0.9,''',letras{t},''',''Units'',''normalized'');'])
            
            axes('Position',[left bottom2 width height2])
            eval(['plot(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'*100 - MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'*100,''k'',''linewidth'',1.5);'])
            hold on
            plot([-88 88],[0 0],'k--')
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            set(gca,'YAxisLocation','right');
            eval(['axis([-88 88 ',axish2min{v},' ',axish2max{v},'])'])
            ylabel('Difference','fontsize',14)
            xlabel('Latitude','fontsize',14)
            
            tightfig
            
        else
            axes('Position',[left bottom1 width height1]);
            eval(['H(1) = shadedErrorBar(lat,MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},'_std, {@mean, @(x) 1*std(x)}, {''-b'', ''LineWidth'', 1.5}, 0);'])
            hold on
            eval(['H(2) = shadedErrorBar(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},'_std, {@mean, @(x) 1*std(x)}, {''-r'', ''LineWidth'', 1.5}, 0);'])
            plot([-88 88],[0 0],'k--')
            if t == 3;
                legend([H(1).mainLine, H(2).mainLine],'<CONTROL>','<TSI10p>','location','best');
            end
            %             pbaspect([1 .3 1])
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            eval(['axis([-88 88 ',axish1min{v},' ',axish1max{v},'])'])
%             eval(['title(''',ttl{v},' - ',trimestre{t},''',''fontsize'',11)'])
            eval(['h_label=ylabel({''',ylbl2{v},''';''',ylbl3{v},'''},''fontsize'',14,''visible'',''on'');'])
            set(gca, 'XTickLabel', [],'XTick',[])
            
            if v == 3;
                eval(['set(gca,''ytick'',[',axish1min{v},'+25:25:',axish1max{v},']);'])
            elseif v == 5;
                eval(['set(gca,''ytick'',[',axish1min{v},'+50:50:',axish1max{v},']);'])
            elseif v == 6;
                eval(['set(gca,''ytick'',[',axish1min{v},'+10:10:',axish1max{v},']);'])
            elseif v == 7;
                eval(['set(gca,''ytick'',[',axish1min{v},'+5:5:',axish1max{v},']);'])
            elseif v == 8;
                eval(['set(gca,''ytick'',[',axish1min{v},'+2:2:',axish1max{v},']);'])
            elseif v == 9;
                eval(['set(gca,''ytick'',[',axish1min{v},'+0.5:0.5:',axish1max{v},']);'])
            end
            
            eval(['text(0.05,0.9,''',letras{t},''',''Units'',''normalized'');'])
            
            
            axes('Position',[left bottom2 width height2])
            eval(['plot(lat,MEDIA_BCN_constsol_10p_',cesm2{v},'_',trimestre{t},' - MEDIA_BCN_2_',cesm2{v},'_',trimestre{t},',''k'',''linewidth'',1.5);'])
            hold on
            plot([-88 88],[0 0],'k--')
            set(gca,'XTick',[-100 -90 -60 -30 0 30 60 90 100]);
            set(gca,'XTickLabel',{'100°S' '90°S' '60°S' '30°S' '0°' '30°N' '60°N' '90°N' '100°N'})
            set(gca,'YAxisLocation','right');
            eval(['axis([-88 88 ',axish2min{v},' ',axish2max{v},'])'])
            ylabel('Difference','fontsize',14)
            xlabel('Latitude','fontsize',14)
            
            tightfig
            
        end
        
        
        
        cd ../figuras
        
        eval(['namest = ''cesm_',cesm2{v},'_',trimestre{t},'_bias_junto_LS'';'])
        saveas(gcf,namest,'fig')
        print(namest,'-dpng','-r300');
        close
        
        cd ../dados
        
    end
end
    

