function [BER, trial_cnt]=calc_BER(path_bit_ref, path_bit_rcv)

% % %% BER Calculation
fbit_snd = fopen(path_bit_ref, 'rt');
fbit_rcv = fopen(path_bit_rcv, 'rt');

bit_cnt = 0;
err_cnt = 0;
trial_cnt = 0;
while (1)
    [bit_snd,len] = fread(fbit_snd, 1, 'char');
    [bit_rcv] = fread(fbit_rcv, 1, 'char');
    
    if len == 0 
        break;
    else
%         if bit_rcv == 10
%             fseek(fbit_snd, 0, 'bof');
%             trial_cnt = trial_cnt + 1;
%         else
            bit_cnt = bit_cnt + 1;
            if ~isempty(bit_rcv)
                if bit_snd ~= bit_rcv 
                    err_cnt = err_cnt + 1;
                end
            else
                err_cnt = err_cnt + 1;
            end
%         end
    end
end

if bit_cnt == 0
    BER = 1;
else
    BER = (err_cnt / bit_cnt);
end

% if BER > 0.5
%     BER = NaN;
% end
disp(['Number of Trial: ' num2str(trial_cnt)]);
disp(['Bit Error Rate = ' num2str(BER)]);