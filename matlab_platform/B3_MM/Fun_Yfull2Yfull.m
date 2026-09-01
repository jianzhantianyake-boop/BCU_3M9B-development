% =========================================================================
% 教学操作说明：根据故障状态修改完整网络导纳矩阵的函数
% 使用方法：
%  按以下函数签名调用：function [Y_full_modified,Transform]=Fun_Yfull2Yfull(Y_full,pfdata,faultflag)
% 参数：
%  以签名中的输入变量为准。状态向量、导纳矩阵和参数结构体的维度及排序必须与初始化结果一致。
% 返回 / 工作区结果：
%  以签名左侧的输出变量为准；结果通常供上层 CCT、轨迹、能量或稳定区域脚本继续使用。
% 步骤：
%  1. 接收上层脚本提供的状态与参数。 2. 按原始方程或迭代规则完成本步骤。 3. 返回残差、状态、能量、矩阵或收敛标志。
% 单位：
%   角度通常为 rad，角速度为 rad/s，时间为 s，功率、电压和导纳通常为 pu；
%   个别中间变量为无量纲标志、迭代次数或矩阵索引，具体以变量定义为准。
% 前置条件：
%  仅在其上层流程已完成初始化后调用；不要仅凭单次函数输出推断整个系统的稳定性。
% 研究与验证边界：
%   本次只更新教学注释，不改变原始方程、参数、判据、求解器或绘图逻辑；
%   MATLAB 原生运行、收敛性和物理结论仍须在目标 MATLAB 环境中实际核验。
% =========================================================================
%% function: modify structure preserved admittance matrix
function [Y_full_modified,Transform]=Fun_Yfull2Yfull(Y_full,pfdata,faultflag)
%% block matrix
    if(faultflag==0)
        nbus=pfdata.bus.num;
        no_faultbus=nbus+1;
    else
        nbus=pfdata.bus.num-1;
        no_faultbus=faultflag(2);
    end
    ngen=pfdata.bus.numgen;
    Ynn=zeros(ngen,ngen);
    Ynr=zeros(ngen,nbus-ngen);
    Yrn=zeros(nbus-ngen,ngen);
    Yrr=zeros(nbus-ngen,nbus-ngen);    
    % Ynn
        for i=1:ngen      
            for j=1:ngen
                no_genself=pfdata.gen.no(i,1);
                no_geninte=pfdata.gen.no(j,1);
                if(no_genself>no_faultbus)
                    no_genself=no_genself-1;
                end
                if(no_geninte>no_faultbus)
                    no_geninte=no_geninte-1;
                end
                Ynn(i,j)=Y_full(no_genself,no_geninte);
            end
        end
    % Ynr
        no_withoutgen=zeros(nbus-ngen,1);
        k=1;
        for i=1:nbus
            flag_gen=0;
            for j=1:ngen
                genno = pfdata.gen.no(j);
                if(genno>no_faultbus)
                    genno = genno-1;
                end
                if(i==genno)
                    flag_gen=1;
                end
            end
            if(flag_gen==0)
                no_withoutgen(k,1)=i;
                k=k+1;
            end
            flag_gen=0; 
        end
        for i=1:ngen      
            for j=1:(nbus-ngen)
                no_genself=pfdata.gen.no(i,1);
                if(no_genself>no_faultbus)
                    no_genself=no_genself-1;
                end
                no_loadinte=no_withoutgen(j,1);
                Ynr(i,j)=Y_full(no_genself,no_loadinte);
            end
        end    
    % Yrn
        for i=1:(nbus-ngen)      
        for j=1:ngen
            no_loadself=no_withoutgen(i,1);
            no_geninte=pfdata.gen.no(j,1);
            if(no_geninte>no_faultbus)
                no_geninte=no_geninte-1;
            end
            Yrn(i,j)=Y_full(no_loadself,no_geninte);
        end
        end   
    % Yrr
        for i=1:(nbus-ngen)      
        for j=1:(nbus-ngen) 
            no_loadself=no_withoutgen(i,1);
            no_loadinte=no_withoutgen(j,1);
            Yrr(i,j)=Y_full(no_loadself,no_loadinte);
        end
        end    
%% modified matrix
    Y_full_modified=[Ynn,Ynr;Yrn,Yrr];
    k=1;
     for i=1:pfdata.bus.num
            flag_gen=0;
            for j=1:ngen
                genno = pfdata.gen.no(j);
                if(i==genno)
                    flag_gen=1;
                end
            end
            if(flag_gen==0 && i~=no_faultbus)
                no_withoutgen(k,1)=i;
                k=k+1;
            end
    end
    Transform = [pfdata.gen.no;no_withoutgen];
end
