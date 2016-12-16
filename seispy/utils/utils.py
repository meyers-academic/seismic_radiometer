def set_channel_vector(letter):
    if letter.upper() == 'Z':
        return [0, 0, 1]
    if letter.upper() == 'N':
        return [0, 1, 0]
    if letter.upper() == 'E':
        return [1, 0 ,0]
