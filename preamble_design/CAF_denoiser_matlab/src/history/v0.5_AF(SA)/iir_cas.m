function [out, buff_out] = iir_cas(in,SOS,G, f_len,buff_in)
    out = zeros(f_len,1);
    cas_len = size(SOS,1);
    buff = buff_in;

    for i = 1 : f_len
          cas_in = in(i) * G(1);
        for j = 1 : cas_len
     %       cas_out = filter(SOS(j,1:3),SOS(j,4:6),cas_in);

             a_sum = cas_in - SOS(j,5)*buff(j,1) - SOS(j,6)*buff(j,2);
             cas_out = (a_sum + SOS(j,2)*buff(j,1) + SOS(j,3)*buff(j,2));
             cas_in = cas_out * G(1+j);
             buff(j,2) = buff(j,1);
             buff(j,1) = a_sum;
        end
        out(i) = cas_in;
    end
    buff_out = buff;
end