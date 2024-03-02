#! /bin/bash
echo "This is the easy-understanding interacting script to execute run.py."
echo "You can also directly run the Python script"
read -e -p "Please input path of gfile:" gfilepath
read -e -p "Please specify the output spdata file path (default: ./spdata.dat): " outputpath
outputpath=${outputpath:-./spdata.dat}
read -e -p "Please input npsi (=lsp in GTC and spdata, default: 91):" lsp
lsp=${lsp:-91}
read -e -p "Please input ntheta (=lst in GTC, =lst+1 in spdata, default: 122):"
lst=${lst:-122}
read -e -p "Please input the ratio of last mapping surface to the last closed flux surface (default: 0.995): " psimax
psimax=${psimax:-0.995}

while true
do
    read -e -p "Do you want to save the consistency checking figures (y/n/none, default: n. none for no figure)?" savefig
    savefig=${savefig:-n}
    
    if [ ${savefig} = "n" ]
    then
        figs="show"
        break
    elif [ ${savefig} = "y" ]
    then
        figs="save"
        break
    elif [ ${savefig} = "none" ]
    then
        figs="none"
        break
    else
        echo "only y/n/none are accepted"
    fi
done

cmd=(python ./run.py
    --input=${gfilepath}
    --output=${outputpath}
    --lsp=${lsp}
    --lst=${lst}
    --psimax=${psimax}
    --figs=${figs}
)

echo ${cmd[@]}
"${cmd[@]}"