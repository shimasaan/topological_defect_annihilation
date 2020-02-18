% #######################################################################################################
% ���]�ڎ����@+-1�g�|���W�J�����ׂ̑Ώ��ŃV�~�����[�V����
% 2019/12/10
% ���R�[�h���璷�Ȃ̂͋C�ɂ��Ă͂����Ȃ�
% #######################################################################################################

% �ϐ��錾
L = 50; %�V�X�e���T�C�Y
a = 0.10; %x�����ւ̃o�C�A�X�W��
u = zeros(L,L); %�z����(x����)
v = zeros(L,L); %�z����(y����)
x = 1:L; %�T�C�g�̈ʒu(x����)
y = (1:L).'; %�T�C�g�̈ʒu(y����)
for i = 1:L-1
    x = [x;1:L];
    y = [y (1:L).'];
end

%% �g�|���W�J�����ׂ���������Ԃ̐���
Xplus = 1*L/4 + L/7.6; %(+1���ׂ�x���W�By���W�͒��S�ɌŒ�)
Xminus = 3*L/4 - L/7.6; %(-1���ׂ�x���W)
for i = 1:L
    for j = 1:L %(+1��)
        u(i,j) = (j-Xplus)/sqrt((i-L/2)^2+(j-Xplus)^2); %�璷���������ǋC�ɂ��Ȃ�
        v(i,j) = (i-L/2)/sqrt((i-L/2)^2+(j-Xplus)^2);
        u(i,j) = (u(i,j)-sign(i-L/2)*(i-L/2)*a); %�z����S�̂ɐ��������փo�C�A�X��������
        z = u(i,j)^2+v(i,j)^2; %�K�i��
        u(i,j) = u(i,j)/sqrt(z);
        v(i,j) = v(i,j)/sqrt(z);
    end
    for j = L/2:L %(-1��)
        v(i,j) = (i-L/2)/sqrt((i-L/2)^2+(j-Xminus)^2);
        u(i,j) = -(j-Xminus)/sqrt((i-L/2)^2+(j-Xminus)^2);
        u(i,j) = (u(i,j)-sign(i-L/2)*(i-L/2)*a);
        z = u(i,j)^2+v(i,j)^2;
        u(i,j) = u(i,j)/sqrt(z);
        v(i,j) = v(i,j)/sqrt(z);
    end
end
% ���E�ł͊��S�ɐ����Ȕz����ɂ���
u(1,:) = -ones(1,L);
v(1,:) = zeros(1,L);
u(:,1) = -ones(L,1);
v(:,1) = zeros(L,1);
u(L,:) = -ones(1,L);
v(L,:) = zeros(1,L);
u(:,L) = -ones(L,1);
v(:,L) = zeros(L,1);
% quiver(x,y,u,v,'ShowArrowHead','on','LineWidth',1.5)

%% ����ۑ��`���̐錾
mov = VideoWriter('mov.mp4','MPEG-4');
mov.FrameRate = 10;
open(mov);

%���������̕\��
quiver(x,y,u,v,'ShowArrowHead','off','LineWidth',1.5,'AutoScaleFactor',0.7);
xlim([0 L]);
ylim([0 L]);
xlabel(['t = ',num2str(0)],'FontSize',12); %�O���t�^�C�g���Ƃ��āA���ԃ��x����\��
img = getframe(gcf);
writeVideo(mov, img);

%�ɘa�ߒ�
T = 250;
for t = 1:T
    for t2 = 1:L*L*2.5

        %�A�b�v�f�[�g����T�C�g�������_���ɑI���i�Œ苫�E�j
%         i = 1+randi(L-2);
%         j = 1+randi(L-2);
        %�A�b�v�f�[�g����T�C�g�������_���ɑI���i�������E�j
        i = randi(L);
        j = randi(L);
        
        %�ŋߐڊi�q�̑I���Ɣz����ǂݍ���
        nnu = [u(rem(i,L)+1,j) u(rem(i+L-2,L)+1,j) u(i,rem(j,L)+1) u(i,rem(j+L-2,L)+1)];
        nnv = [v(rem(i,L)+1,j) v(rem(i+L-2,L)+1,j) v(i,rem(j,L)+1) v(i,rem(j+L-2,L)+1)];
        
        %nematic���ݍ�p�̏ꍇ�͎��s(polar�̏ꍇ�͕s�v)
%         for k = 1:4
%            nnu(1,k) = sign(u(i,j)*nnu(1,k)+v(i,j)*nnv(1,k))*nnu(1,k);
%            nnv(1,k) = sign(u(i,j)*nnu(1,k)+v(i,j)*nnv(1,k))*nnv(1,k);
%         end
        
        %���ω������z���׃N�g���̍X�V
        u(i,j) = mean(nnu);
        v(i,j) = mean(nnv);
        z = sqrt(u(i,j)^2+v(i,j)^2);
        u(i,j) = u(i,j)/z;
        v(i,j) = v(i,j)/z;     
    end
    
    %�t���[����\��
    quiver(x,y,u,v,'ShowArrowHead','off','LineWidth',1.5,'AutoScaleFactor',0.7);
    xlim([0 L]);
    ylim([0 L]);
    xlabel(['t = ',num2str(t)],'FontSize',12); %�O���t�^�C�g���Ƃ��āA���ԃ��x����\��
    img = getframe(gcf);
    writeVideo(mov, img);
end
close(mov)