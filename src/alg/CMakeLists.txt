aux_source_directory(. DIR_ALG_SRCS)
add_library(alg ${DIR_ALG_SRCS})

# 指定 MathFunctions 库的安装路径
install (TARGETS alg DESTINATION ./bin)