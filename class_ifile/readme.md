# Descrição
A classe *class_ifile* foi desenvolvida para simplificar a leitura de parâmetros de entrada em códigos que usam a linguagem Fortran 2008.

# Como usar

## Prepare o arquivo de entrada

O arquivo de entrada deve seguir o seguinte padrão

```
valor FS nome_da_variável FS descrição da variável
```

onde FS é o delimitador de campo. Qualquer linha que não siga este padrão é tratada como comentário.

## Inclua o módulo
Use o comando
```
use mod_class_ifile
```
para incluir o módulo da classe.

## Instancie e inicialize o objeto

Crie uma instância da classe, por exemplo:
```
type(class_ifile) :: ifile
```

Inizialize o objeto:
```
call ifile%init("name of the input file","field delimiter")
```

Carregue o arquivo de entrada no objeto:
```
call ifile%load()
```

## Leia os parâmetros de interesse

Para obter o valor de uma variável *VAR* cujo nome é *VARNAME* no arquivo de entrada, faça
```
call ifile%get_value(VAR,"VARNAME")
```

## Exemplo

Para um exemplo concreto de uso da classe, analise o arquivo *example.f90*. 

Para compilar o exemplo usando o compilador *gfortran*, execute a seguinte linha em um terminal:
```
gfortran class_ifile.f90 example.f90 -o example.exe
```

Para executar o programa, execute no terminal o seguinte comando:
```
./example.exe
```