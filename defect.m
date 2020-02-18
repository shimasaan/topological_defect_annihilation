% #######################################################################################################
% 相転移実験　+-1トポロジカル欠陥の対消滅シミュレーション
% 2019/12/10
% ※コードが冗長なのは気にしてはいけない
% #######################################################################################################

% 変数宣言
L = 50; %システムサイズ
a = 0.10; %x方向へのバイアス係数
u = zeros(L,L); %配向場(x成分)
v = zeros(L,L); %配向場(y成分)
x = 1:L; %サイトの位置(x成分)
y = (1:L).'; %サイトの位置(y成分)
for i = 1:L-1
    x = [x;1:L];
    y = [y (1:L).'];
end

%% トポロジカル欠陥を持つ初期状態の生成
Xplus = 1*L/4 + L/7.6; %(+1欠陥のx座標。y座標は中心に固定)
Xminus = 3*L/4 - L/7.6; %(-1欠陥のx座標)
for i = 1:L
    for j = 1:L %(+1側)
        u(i,j) = (j-Xplus)/sqrt((i-L/2)^2+(j-Xplus)^2); %冗長かもだけど気にしない
        v(i,j) = (i-L/2)/sqrt((i-L/2)^2+(j-Xplus)^2);
        u(i,j) = (u(i,j)-sign(i-L/2)*(i-L/2)*a); %配向場全体に水平方向へバイアスをかける
        z = u(i,j)^2+v(i,j)^2; %規格化
        u(i,j) = u(i,j)/sqrt(z);
        v(i,j) = v(i,j)/sqrt(z);
    end
    for j = L/2:L %(-1側)
        v(i,j) = (i-L/2)/sqrt((i-L/2)^2+(j-Xminus)^2);
        u(i,j) = -(j-Xminus)/sqrt((i-L/2)^2+(j-Xminus)^2);
        u(i,j) = (u(i,j)-sign(i-L/2)*(i-L/2)*a);
        z = u(i,j)^2+v(i,j)^2;
        u(i,j) = u(i,j)/sqrt(z);
        v(i,j) = v(i,j)/sqrt(z);
    end
end
% 境界では完全に水平な配向場にする
u(1,:) = -ones(1,L);
v(1,:) = zeros(1,L);
u(:,1) = -ones(L,1);
v(:,1) = zeros(L,1);
u(L,:) = -ones(1,L);
v(L,:) = zeros(1,L);
u(:,L) = -ones(L,1);
v(:,L) = zeros(L,1);
% quiver(x,y,u,v,'ShowArrowHead','on','LineWidth',1.5)

%% 動画保存形式の宣言
mov = VideoWriter('mov.mp4','MPEG-4');
mov.FrameRate = 10;
open(mov);

%初期条件の表示
quiver(x,y,u,v,'ShowArrowHead','off','LineWidth',1.5,'AutoScaleFactor',0.7);
xlim([0 L]);
ylim([0 L]);
xlabel(['t = ',num2str(0)],'FontSize',12); %グラフタイトルとして、時間ラベルを表示
img = getframe(gcf);
writeVideo(mov, img);

%緩和過程
T = 250;
for t = 1:T
    for t2 = 1:L*L*2.5

        %アップデートするサイトをランダムに選択（固定境界）
%         i = 1+randi(L-2);
%         j = 1+randi(L-2);
        %アップデートするサイトをランダムに選択（周期境界）
        i = randi(L);
        j = randi(L);
        
        %最近接格子の選択と配向場読み込み
        nnu = [u(rem(i,L)+1,j) u(rem(i+L-2,L)+1,j) u(i,rem(j,L)+1) u(i,rem(j+L-2,L)+1)];
        nnv = [v(rem(i,L)+1,j) v(rem(i+L-2,L)+1,j) v(i,rem(j,L)+1) v(i,rem(j+L-2,L)+1)];
        
        %nematic相互作用の場合は実行(polarの場合は不要)
%         for k = 1:4
%            nnu(1,k) = sign(u(i,j)*nnu(1,k)+v(i,j)*nnv(1,k))*nnu(1,k);
%            nnv(1,k) = sign(u(i,j)*nnu(1,k)+v(i,j)*nnv(1,k))*nnv(1,k);
%         end
        
        %平均化した配向べクトルの更新
        u(i,j) = mean(nnu);
        v(i,j) = mean(nnv);
        z = sqrt(u(i,j)^2+v(i,j)^2);
        u(i,j) = u(i,j)/z;
        v(i,j) = v(i,j)/z;     
    end
    
    %フレームを表示
    quiver(x,y,u,v,'ShowArrowHead','off','LineWidth',1.5,'AutoScaleFactor',0.7);
    xlim([0 L]);
    ylim([0 L]);
    xlabel(['t = ',num2str(t)],'FontSize',12); %グラフタイトルとして、時間ラベルを表示
    img = getframe(gcf);
    writeVideo(mov, img);
end
close(mov)