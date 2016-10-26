function s = sample_cut(s,len)

if length(s) > len
    s = s(1:len);
end

if length(s) < len
    while 1
        s = [s;s];
        if length(s) >= len
            s = s(1:len);
            break;
        end
    end
end

end