rm log.log
for i in */
do
  cd $i
  if [ -a doney.log ]; then
    cd ..
    echo $i
    echo "already named moving on"
    continue
  fi
  mkdir -p dicom
  mv *.dcm dicom
  cd dicom
  #var=$(find *.dcm -print -quit)
  for var in *.dcm
  do
    gdcmraw -i $var -t 20,1002 -o slicenumber.log
    slices=$(cat slicenumber.log)
    #slices="2728"
    #use 2000 as number of slices to check the quality of DWI acquisitions
    #if [ $slices -gt "2000" ]; then
    if [ $slices -gt "100" ]; then
      echo $i
      echo $var
      echo "###############################" >> ../../log.log
      echo $var >> ../../log.log
      gdcmdump $var | grep "Patient ID" >> ../../log.log
      gdcmraw -i $var -t 10,20 -o ../../ID.log
      ID=$(cat ../../ID.log)
      gdcmdump $var | grep "Modality" >> ../../log.log
      gdcmdump $var | grep "Acquisition Date" >> ../../log.log
      gdcmraw -i $var -t 8,22 -o ../../Date.log
      Date=$(cat ../../Date.log)
      gdcmdump $var | grep "Images in Acquisition" >> ../../log.log
      gdcmdump $var | grep "Last image number used" >> ../../log.log
      #cd ..
      outputname=${ID% TWH }
      outputname+="-"
      outputname+=$Date
      #mv $i $outputname
      break
    fi
  done

  cd ..
  cd ..
  if [ -n "$outputname" ]; then
    mv $i $outputname
  fi
  rm Date.log
  rm ID.log
done

for dirs in *-*/
do
  cd $dirs
  echo $dirs
  if [ -a doney.log ]; then
    cd ..
    #echo $dirs
    echo "already done moving on."
    continue
  fi
  cd dicom
  for file in *.dcm; do gdcmraw -i $file -t 8,103e -o name.log; echo $(cat name.log) >> names.txt; done
  cat names.txt | sort | uniq >> new.txt
  while read filetype
  do
    mkdir -p ../"$filetype"
  done < new.txt
  for file in *.dcm
  do
    gdcmraw -i $file -t 8,103e -o name.log
    name=$(cat name.log)
    while read filetype
    do
      if [[ $name = "$filetype " ]]; then
        ln -s -t ../"$filetype" ../dicom/$file
      elif [[ $name = "$filetype" ]]; then
        ln -s -t ../"$filetype" ../dicom/$file
      fi
    done < new.txt
  done
  names.txt name.log slicenumber.log
  cd ..
  while read filetype
  do
    cd "$filetype"
    for i in *.dcm
    do
      gdcmconv $i $i -w
    done
    cd ..
    dcm2nii -d N -e N -f Y -i N -p N $filetype
    rm "$filetype"/*o*.nii.gz
    #rmdir "$filetype"
  done < dicom/new.txt
  rm dicom/new.txt
  echo $dirs >> doney.log
  cd ..
  echo "##############"
done
