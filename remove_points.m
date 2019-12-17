function [Snew,Prm_new,Pmm_new,lmkinfo,feats_list]=remove_points(Snew,Prm_new,Pmm_new,lmkinfo,feats_list)

tbr=find(lmkinfo.counter_meas./lmkinfo.counter_prop<0.5 & lmkinfo.counter_prop>5);

feats_list(tbr,:)=[];
% if ~isempty(tbr)
%     disp('wait');
% end
lmkinfo.counter_meas(tbr)=[];
lmkinfo.counter_prop(tbr)=[];


rows=reshape(transpose((tbr-1)*3+[1:3]),[],1);

Snew(tbr,:)=[];
Prm_new(:,rows)=[];
Pmm_new(rows,:)=[];
Pmm_new(:,rows)=[];






end