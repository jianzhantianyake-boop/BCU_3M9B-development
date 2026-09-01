% =========================================================================
% 教学操作说明：由线路 RXB 数据构造完整导纳矩阵的函数
% 使用方法：
%  按以下函数签名调用：function Yfull=Fun_RXB2Yfull(RXB,pfdata)
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
%% function: transfer RXB (from matpower results) to admittance matrix (full order)
function Yfull=Fun_RXB2Yfull(RXB,pfdata)
    Yii0=zeros(pfdata.bus.num,2);
    Yij_re=zeros(pfdata.bus.num,pfdata.bus.num);
    Yij_im=zeros(pfdata.bus.num,pfdata.bus.num);
    for i=1:pfdata.bus.num  % search for the bus i
        for k=1:size(RXB,1)  % scan the k line
            if(RXB(k,1)==i||RXB(k,2)==i)
                Yii0(i,1)=Yii0(i,1)+RXB(k,3)/(RXB(k,3)^2+RXB(k,4)^2);
                Yii0(i,2)=Yii0(i,2)+RXB(k,5)/2-RXB(k,4)/(RXB(k,3)^2+RXB(k,4)^2);
                if(RXB(k,1)==i)
                    j=RXB(k,2);
                else
                    j=RXB(k,1);
                end
                Yij_re(i,j)=-1*RXB(k,3)/(RXB(k,3)^2+RXB(k,4)^2);
                Yij_im(i,j)=RXB(k,4)/(RXB(k,3)^2+RXB(k,4)^2);
                Yij_re(j,i)=Yij_re(i,j);
                Yij_im(j,i)=Yij_im(i,j);
            end
        end
    end
    Yij=complex(Yij_re,Yij_im);        
    Yii=complex(Yii0(:,1),Yii0(:,2));
    Yii_mat=diag(Yii);
    Yfull=Yij+Yii_mat;
end
