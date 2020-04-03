import pandas as pd
import json
import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped


class Patient:
    """ Класс для хранения данных пациента """

    def __init__(self, pid, data_path='data/hivevo', reference_path='data/references', verbose=False):
        """
        Args:
            pid (str): ID пациента
            data_path (str): путь до данных
            verbose (bool): говорливость кода
            reference_path (str): путь, куда будет произведена загрузка референсов
        """
        self.id = pid
        self.reference = Reference(reference_path, patient=pid)
        self.status = None
        self.data_path = data_path
        self.regions = None
        self.verbose = verbose

        # сразу записываем регионы
        self.init_regions()

    def init_regions(self):
        """
        Инициализирует данные о регионах
        """
        regions = pd.DataFrame()

        # цикл по возможным регионам
        for name in ["V3", "PR", "psi", "vpr", "vpu", "p1", "p2", "p6", "p7", "p15", "p17", "RRE"]:
            filename = f'hivevo_{self.id}_{name}.json'

            with open(os.path.join(self.data_path, filename)) as f:
                region_file = json.load(f)

            # если не смогли скачать данные, то пишем об этом
            if type(region_file) == dict:
                if self.verbose:
                    print(f'No haplotype for patient {self.id} for region {name}')
                continue

            # временный датасет
            region = pd.DataFrame(data=region_file)
            region.columns = ["days", "frequency", "sequence", "name", "description"]

            # вырезаем число прочтений
            region['nreads'] = region.description.str.split(': ').str[-1]

            # заменяем неинформативную колонку на название региона
            region['name'] = name

            region['translated'] = region.sequence.apply(lambda x: Seq(x,
                                                                       Gapped(IUPAC.unambiguous_dna)
                                                                       ).ungap().translate())

            # удаляем ненужную
            region.drop(['description'], axis=1, inplace=True)
            regions = pd.concat([regions, region], ignore_index=True)

        # теперь имеем объект регионов!
        # сюда дописываем данные по референсам, они под днем 0
        self.regions = Region(pd.concat([regions, self.reference.reference_df.drop(['id'], axis=1)], ignore_index=True))


class Reference:
    """ Класс для работы с данными референсов """

    def __init__(self, reference_path, patient=None):

        self.reference_df = pd.DataFrame()

        # сначала получаем все названия файлов
        reference_list = [x for x in os.listdir(reference_path) if 'reference' in x]

        # если нужно было выбрать одного пациента, то оставляем только соответствующие названия
        if patient:
            reference_list = [x for x in reference_list if f'_{patient}.' in x]

        # собираем данные из json-ов
        for name in reference_list:
            with open(os.path.join(reference_path, name)) as f:
                json_file = json.load(f)

                # удаляем ненужные колонки
                t = pd.DataFrame(data=json_file).drop(['name', 'description'], axis=1)

                # делим объединённые колонки на отдельные
                t = pd.concat([t.drop(['features'], axis=1), t.features.apply(pd.Series)], axis=1)

                # переводим в простой список
                t.location = t.location.apply(pd.Series)

                # вырезаем последовательность
                t['region_seq'] = t.apply(lambda x: x.seq[x.location[0]:x.location[1]].strip(), axis=1)

                # переименовываем для единообразия
                t.rename(mapper={'region_seq': 'sequence', 'seq': 'full_reference'}, axis=1, inplace=True)
                t['translated'] = t.sequence.apply(lambda x: Seq(x, Gapped(IUPAC.unambiguous_dna)).ungap().translate())
                self.reference_df = pd.concat([self.reference_df, t], ignore_index=True)

        # оставляем только нужные колонки
        self.reference_df = self.reference_df[['sequence', 'name', 'translated', 'id']]

        # остальное дописываем вручную (это будет использоваться в дереве)
        self.reference_df['days'] = 0
        self.reference_df['frequency'] = 100
        self.reference_df['nreads'] = 1

        # если это делалось для одного пациента, то сразу отдаём объект с его регионами
        if patient:
            self.region = Region(self.reference_df)

    def get_patient(self, patient, region=None):
        """
        Функция для извлечения данных о пациенте или пациенте и регионе, если используем класс отдельно
        Args:
            patient (str): ID пациента
            region (str): регион

        Returns:
            pd.DataFrame: отфильтрованные данные
        """
        if region:
            return self.reference_df.loc[(self.reference_df.id == f'reference_{patient}') &
                                         (self.reference_df.name == region)]
        else:
            return self.reference_df.loc[self.reference_df.id == f'reference_{patient}']


class Region:
    """ Класс для взаимодействия с регионами """

    def __init__(self, data):
        """
        Args:
            data (pd.DataFrame): данные пациента
        """
        self.region = data
        self.keys = data.name.unique()

        # TODO: скоринг сразу
        self.score()

    def __getitem__(self, key):
        return self.region.loc[self.region.name == key]

    def __len__(self):
        return len(self.keys)

    def __repr__(self):
        return repr(self.region)

    def score(self):
        """
        Скоринг по Кидера факторам и свойствам аминокислот
        """
        pass

    def major(self, region=None):
        """
        Получение наиболее частых гаплотипов за день

        Args:
            region (str): регион

        Returns:
            pd.DataFrame: отфильтрованные данные
        """

        # если используем с данными пациента (не референс)
        if 'days' in self.region.columns:

            # если нужен отдельный регион
            if region:
                return self.region.loc[self.region.name ==
                                       region].sort_values(['frequency'], ascending=[False]).groupby(['days']).first()

            return self.region.sort_values(['frequency'], ascending=[False]).groupby(['name', 'days']).first()

        # если данные о референсе
        return None

