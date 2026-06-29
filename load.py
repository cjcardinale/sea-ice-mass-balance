#make a thickness threshold variable at the beginning (maybe set it to 10)
#e.g., thick_thresh=10 #m
import xarray as xr
import cf_xarray as cfxr
import xesmf as xe
import intake
from intake_esm import DerivedVariableRegistry
import regionmask
import numpy as np
import pandas as pd
from dask import delayed,compute
from tqdm.auto import tqdm
import geopandas as gp
from xmip.preprocessing import rename_cmip6, promote_empty_dims, broadcast_lonlat, replace_x_y_nominal_lat_lon, combined_preprocessing, correct_lon, parse_lon_lat_bounds, correct_coordinates, maybe_convert_bounds_to_vertex, maybe_convert_vertex_to_bounds, sort_vertex_order, fix_metadata, correct_units
from importlib.resources import files
import os
import warnings
from pyesgf.search import SearchConnection
warnings.filterwarnings('ignore')
from io import StringIO  
import sys
from copy import deepcopy
from pyproj import Geod
from scipy.interpolate import griddata
from dask.diagnostics import ProgressBar
import intake_esgf
import glob
#intake_esgf.conf.set(all_indices=True) 
#intake_esgf.conf.set(all_indices=True) 
#intake_esgf.conf.set(indices={"esgf.nci.org.au":False})
#intake_esgf.conf.set(indices={"esgf.nci.org.au":False})
#intake_esgf.conf.set(indices={"esg-dn1.nsc.liu.se":False})
#intake_esgf.conf.set(indices={"esgf.ceda.ac.uk":False})
#intake_esgf.conf.set(indices={"esgf-data.dkrz.de":False})
#intake_esgf.conf.set(indices={"esgf-node.ipsl.upmc.fr":False})
#intake_esgf.conf.set(indices={"esgf-node.ornl.gov":True})
#intake_esgf.conf.set(indices={"esgf-node.llnl.gov":False})
#intake_esgf.conf.set(indices={"ESGF2-US-1.5-Catalog":True})
#intake_esgf.conf.set(indices={"anl-dev":True})
#intake_esgf.conf.set(indices={"ornl-dev":True})
#intake_esgf.conf["break_on_error"] = False
intake_esgf.conf.set(all_indices=True, break_on_error=False)
from intake_esgf import ESGFCatalog
from intake_esgf.exceptions import NoSearchResults
from globus_sdk.services.search.errors import SearchAPIError
from collections import ChainMap, defaultdict  
import cftime
from urllib.parse import urlparse
import requests
import tempfile
import time
import random
import xml.etree.ElementTree as ET

# THREDDS OPeNDAP endpoint for CESM2-LE ice monthly time-series files
_CESM2LE_THREDDS_BASE = "https://tds.gdex.ucar.edu/thredds"
_CESM2LE_THREDDS_PATH = "/files/d651056/CESM2-LE/ice/proc/tseries/month_1"

# Session-level cache for THREDDS HTML catalog pages (keyed by URL)
_THREDDS_CATALOG_CACHE: dict = {}

def _griddata(arr, xi, method: str):
    ar1d = arr.ravel()
    valid = np.isfinite(ar1d)
    if valid.all():
        return arr
    return griddata(
        points=tuple(x[valid] for x in xi),
        values=ar1d[valid],
        xi=xi,
        method=method,
        fill_value=np.nan,
    ).reshape(arr.shape)
def interpolate_na(da, dim, method="nearest", use_coordinates=True, keep_attrs=True):
    # Create points only once.
    if use_coordinates:
        coords = [da.coords[d] for d in dim]
    else:
        coords = [np.arange(da.sizes[d]) for d in dim]

    xi = tuple(x.ravel() for x in np.meshgrid(*coords, indexing="ij"))
    arr = xr.apply_ufunc(
        _griddata,
        da,
        input_core_dims=[dim],
        output_core_dims=[dim],
        #output_dtypes=[da.dtype],
        dask="parallelized",
        vectorize=False,
        keep_attrs=keep_attrs,
        kwargs={"xi": xi, "method": method},
    ).transpose(*da.dims)
    return arr

dvr = DerivedVariableRegistry()
@dvr.register(variable='sifb_d', query={'variable_id': ['sithick','sisnthick','siconc']})
def calc_freeboard(ds):
    rho_w = 1026 #kg/m3
    rho_i = 916 #kg/m3
    rho_sn = 330 #kg/m3
    if 'source_id' in ds.attrs:
        if ds.attrs['source_id']=='TaiESM1':
            ds['sithick']=ds['sithick']/ds['siconc']/(ds['siconc']/100)
    H_i = ds.sithick.where(lambda x:np.abs(x)<10)
    H_sn = ds.sisnthick.where(lambda x:np.abs(x)<10)
    ds['sifb_d'] = H_i * ((rho_w - rho_i)/rho_w) - H_sn * (rho_sn/rho_w)
    return ds

@dvr.register(variable='sifb_d2', query={'variable_id': ['sivol','sisnthick','siconc']})
def calc_freeboard2(ds):
    rho_w = 1026 #kg/m3
    rho_i = 916 #kg/m3
    rho_sn = 330 #kg/m3
    H_i = (ds.sivol/(ds.siconc.where(ds.siconc>0)/100)).where(lambda x:np.abs(x)<10)
    H_sn = ds.sisnthick.where(lambda x:np.abs(x)<10)
    ds['sifb_d2'] = H_i * ((rho_w - rho_i)/rho_w) - H_sn * (rho_sn/rho_w)
    return ds

@dvr.register(variable='sifb_d3', query={'variable_id': ['sithick','sisnthick']})
def calc_freeboard(ds):
    rho_w = 1026 #kg/m3
    rho_i = 916 #kg/m3
    rho_sn = 330 #kg/m3
    H_i = ds.sithick.where(lambda x:np.abs(x)<10)
    H_sn = ds.sisnthick.where(lambda x:np.abs(x)<10)
    ds['sifb_d3'] = H_i * ((rho_w - rho_i)/rho_w) - H_sn * (rho_sn/rho_w)
    return ds

@dvr.register(variable='rhoi', query={'variable_id': ['sithick','sisnthick','sifb','siconc']})
def calc_rhoi(ds):
    rho_w = 1026 #kg/m3
    rho_sn = 330 #kg/m3
    if 'source_id' in ds.attrs:
        if ds.attrs['source_id']=='TaiESM1':
            ds['sithick']=ds['sithick']/ds['siconc']/(ds['siconc']/100)
    rho_ice = rho_w - ((rho_w*ds.sifb.where(lambda x:np.logical_and(x!=0, np.abs(x)<10))
                        + rho_sn*ds.sisnthick.where(lambda x:np.logical_and(x!=0, np.abs(x)<10))) / ds.sithick.where(lambda x:np.logical_and(x!=0, np.abs(x)<10)))
    ds['rhoi'] = rho_ice
    return ds

@dvr.register(variable='rhoi2', query={'variable_id': ['sivol','simass','siconc']})
def calc_rhoi2(ds):
    ds['rhoi2'] = ds['simass']/ds['sivol']
    return ds

@dvr.register(variable='sit_d', query={'variable_id': ['sivol','siconc']})
def calc_sit(ds):
    ds['sit_d'] = ds.sivol/(ds.siconc.where(ds.siconc>0)/100)
    return ds

def data_path(file):
    return files('files').joinpath(file)

grid_CESM2 = xr.open_dataset(data_path('CESM2_grid.nc'))
grid_ATL20_nh = xr.open_dataset(data_path('IS2_grid.nc'))
grid_OSISAF_nh = xr.open_dataset(data_path('OSISAF_nh_grid.nc'))
grid_OSISAF_sh = xr.open_dataset(data_path('OSISAF_sh_grid.nc'))
grid_ATL20_sh = xr.open_dataset(data_path('ATL20_grid.nc'))
NH_seaice_regions = gp.read_file(data_path('NSIDC-0780_SeaIceRegions_NH_v1.0.shp'))
SH_seaice_regions =gp.read_file(data_path('NSIDC-0780_SeaIceRegions_SH-NASA_v1.0.shp'))
ATL20_area_NH = xr.open_dataset(data_path('NSIDC0771_CellArea_PS_N25km_v1.0.nc')).cell_area
ATL20_area_SH = xr.open_dataset(data_path('NSIDC0771_CellArea_PS_S25km_v1.0.nc')).cell_area

def add_corners(ds):
        ds2 = ds.copy()
        ds3 = ds.copy()
        if ds.lon.min()>=0:
            ds2['lon'] = ds.lon-180
            ds3['lon'] = ds.lon-180
            lon_correction = 180
        else:
            lon_correction = 0 
        ds2['lon'] = ds2.lon.where(ds2.lon>0)
        ds3['lon'] = ds3.lon.where(ds3.lon<0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds2 = ds2.cf.add_bounds(['lon','lat'])
            ds3 = ds3.cf.add_bounds(['lon','lat'])
        ds3['lon_bounds'] = ds2.lon_bounds.fillna(0) + ds3.lon_bounds.fillna(0) + lon_correction
        ds3['lon'] = ds.lon
        lat_corners = cfxr.bounds_to_vertices(ds3.lat_bounds.chunk(dict(bounds=-1,y=-1,x=-1)), "bounds", order=None)
        lon_corners = cfxr.bounds_to_vertices(ds3.lon_bounds.chunk(dict(bounds=-1,y=-1,x=-1)), "bounds", order=None)
        ds3=ds3.assign_coords(lon_b=lon_corners, lat_b=lat_corners)
        ds3=ds3.rename({'y_vertices':'y_b','x_vertices':'x_b'}).drop_vars(['lat_bounds','lon_bounds'])
        return ds3

def convert_time_old(ds):
    try:
        exp_id = ds.attrs.get("intake_esm_attrs:experiment_id") or ds.attrs.get("experiment_id")
        table_id = ds.attrs.get("table_id")

        # Apply time slicing only if it makes sense
        if exp_id == "historical":
            ds = ds.sel(time=slice("1850-01", "2014-12"))
        elif exp_id == "hist-1950":
            sds = ds.sel(time=slice("1950-01", "2014-12"))
        elif exp_id:
            ds = ds.sel(time=slice("2015-01", "2100-12"))

        # Fix time formatting
        if "time" in ds.coords:
            if table_id == "SIday":
                ds["time"] = pd.to_datetime(ds.time.dt.strftime("%Y-%m-%d").values)
            elif table_id == "SImon":
                ds["time"] = pd.to_datetime(ds.time.dt.strftime("%Y-%m").values)

        # Drop duplicate timestamps
        ds = ds.drop_duplicates("time", keep="first")
        return ds

    except Exception as e:
        #print(f"[convert_time] ⚠️ Error: {e}")
        return ds

def convert_time2_old(ds):
    try:
        exp_id = ds.attrs.get("intake_esm_attrs:experiment_id") or ds.attrs.get("experiment_id")

        if "time" in ds.coords:
            time_vals = ds.time.values

            # ✅ Convert time to datetime64 if needed
            if not np.issubdtype(time_vals.dtype, np.datetime64):
                time_vals = pd.to_datetime(time_vals)

            # ✅ Infer frequency
            inferred_freq = pd.infer_freq(time_vals[:5])  # use first few values

            if inferred_freq is not None:
                if inferred_freq.startswith("M"):  # monthly
                    time_vals = pd.to_datetime(time_vals).to_period("M").to_timestamp()
                elif inferred_freq.startswith("D"):  # daily
                    time_vals = pd.to_datetime(time_vals).to_period("D").to_timestamp()
                else:
                    # fallback: no adjustment
                    time_vals = pd.to_datetime(time_vals)

            ds["time"] = ("time", time_vals)

            # ✅ Time slicing by experiment
            if exp_id == "historical":
                ds = ds.sel(time=slice("1850-01", "2014-12"))
            elif exp_id == "hist-1950":
                ds = ds.sel(time=slice("1950-01", "2014-12"))
            elif exp_id:
                ds = ds.sel(time=slice("2015-01", "2100-12"))

            # ✅ Drop duplicate time entries
            ds = ds.drop_duplicates("time", keep="first")

        return ds

    except Exception as e:
        #print(f"[convert_time] ⚠️ Error: {e}")
        return ds


def sanitize_time(ds):
    """
    Repair time coordinate for CMIP-style datasets.
    Works with cftime and datetime64.
    Handles monthly and daily data.
    """

    ds = ds.copy()

    if "time" not in ds.coords:
        return ds

    t = ds.time.values

    try:
        # If cftime objects, rebuild timestamps
        if isinstance(t[0], cftime.datetime):
            new_time = pd.to_datetime(
                [f"{tt.year:04d}-{tt.month:02d}-{tt.day:02d}" for tt in t]
            )
        else:
            new_time = pd.to_datetime(t)

        ds = ds.assign_coords(time=("time", new_time))

    except Exception as e:
        print(f"[sanitize_time] ⚠️ could not convert time: {e}")
        return ds

    # Sort time
    ds = ds.sortby("time")

    # Remove duplicate timestamps
    _, idx = np.unique(ds.time.values, return_index=True)
    ds = ds.isel(time=np.sort(idx))

    return ds

def convert_time(ds):
    try:
        exp_id = ds.attrs.get("intake_esm_attrs:experiment_id") or ds.attrs.get("experiment_id")

        if "time" in ds.coords:
            ds = sanitize_time(ds)

            if exp_id == "historical":
                ds = ds.sel(time=slice("1850-01-01", "2014-12-31"))
            elif exp_id == "hist-1950":
                ds = ds.sel(time=slice("1950-01-01", "2014-12-31"))
            elif exp_id:
                ds = ds.sel(time=slice("2015-01-01", "2100-12-31"))

        return ds

    except Exception as e:
        print(f"[convert_time] ⚠️ Error: {e}")
        print(ds)
        print(ds.time)
        return ds


def update_member_id(ds):
    if 'variant_label' in ds.variables:
        ds=ds.rename({'variant_label':'member_id'})
    elif 'variant_label' in ds.attrs and 'member_id' not in ds.variables:
        ds=ds.expand_dims({'member_id':[ds.attrs['variant_label']]})
    if 'sub_experiment_id' in ds.variables:
        ds = ds.isel(sub_experiment_id=0).drop_vars('sub_experiment_id')
    ds['member_id'] = (ds.attrs['source_id']+'_') + ds.member_id.astype('object')
    ds['member_id'] = ds.member_id.astype('<U25')
    ds = ds.expand_dims({'experiment_id':[ds.attrs['experiment_id']]})
    ds['experiment_id'] = ds.experiment_id.astype('<U10')
    return ds

def set_chunks(ds,time_chunks):
    return ds.chunk(chunks={'time': time_chunks})

def complete_preprocessing(ds):
    #ds = ds.copy()
    if 'lat' in ds.coords and 'latitude' in ds.coords and 'y' not in ds.coords:
        ds=ds.rename({'lat':'y','lon':'x'}).drop_vars(['lon_bnds','lat_bnds'], errors="ignore")
    ds = rename_cmip6(ds)
    ds = promote_empty_dims(ds)
    ds = correct_coordinates(ds)
    if 'x' in ds.variables or 'lon' in ds.variables:
        ds = broadcast_lonlat(ds)
        ds = correct_lon(ds)
        #ds = parse_lon_lat_bounds(ds)
        ds = sort_vertex_order(ds)
    #ds = correct_units(ds)
    try:
        ds = maybe_convert_bounds_to_vertex(ds)
    except Exception as error:
        pass
    ds = maybe_convert_vertex_to_bounds(ds)
    ds = fix_metadata(ds)
    if 'vertices_latitude' in ds.variables:
        ds = ds.drop_vars(['vertices_latitude','vertices_longitude'],errors='ignore')
    if 'nvertices' in ds.variables:
        ds=ds.rename({'nvertices':'vertex'})
    if 'vertex' in ds.variables and 'y' in ds.variables:
        if ds.vertex.size == 4 and ds.y.size>1:
            if 'lon_verticies' in ds.variables:
                lon_corners = cfxr.bounds_to_vertices(ds.lon_verticies.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
                lat_corners = cfxr.bounds_to_vertices(ds.lat_verticies.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
            if 'lon_bounds' in ds.variables and 'lon_verticies' not in ds.variables:
                lon_corners = cfxr.bounds_to_vertices(ds.lon_bounds.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
                lat_corners = cfxr.bounds_to_vertices(ds.lat_bounds.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
            ds=ds.assign_coords(lon_b=lon_corners, lat_b=lat_corners)
            ds=ds.rename({'y_vertices':'y_b','x_vertices':'x_b'})
            
    return ds

#load data from opendap url, needed to be in a function for list comprehension to work with try/except
def open_ds(file,chunks):
    try: 
        ds = xr.open_dataset(file.opendap_url, chunks={'time': chunks})
        ds = complete_preprocessing(ds)
        return ds
    except Exception as error:
        pass

def calc_areacello_gufunc(lons,lats,lons2x,lats2x,lons2y,lats2y):
    geod = Geod(ellps='WGS84')
    _,_, distEW = geod.inv(lons,lats,lons2x,lats2x)
    _,_, distNS = geod.inv(lons,lats,lons2y,lats2y)
    pixel_area = distEW * distNS
    return pixel_area

def calc_areacello(ds):
    lons=ds.lon
    lons2x=ds.lon.shift(x=1)
    lons2y=ds.lon.shift(y=1)
    lats=ds.lat
    lats2x=ds.lat.shift(x=1)
    lats2y=ds.lat.shift(y=1)
    area = xr.apply_ufunc(
        calc_areacello_gufunc,
        lons,lats,lons2x,lats2x,lons2y,lats2y,
        dask='allowed',
        output_dtypes=[float])
    area = area.fillna(area.min())
    #ds['areacello']=area
    #ds = ds.assign_coords({'areacello':area})
    area['areacello']=area
    return area

def require_all(sub_df,exp,var):
    """
    This function checks if every (`source_id`, `member_id`) combination
    has data for ALL required (`experiment_id`, `variable_id`) pairs.
    """
    if type(exp)==str:
        exp = [exp]
    if type(var)==str:
        var = [var]
    # Check if all required experiment-variable combinations exist for this `source_id`, `member_id`
    for e in exp:
        for v in var:
            if not ((sub_df["experiment_id"] == e) & (sub_df["variable_id"] == v)).any():
                return False  # Missing required data → Remove this `source_id`, `member_id`
    
    return True  # Keep this `source_id`, `member_id`

def filter_missing(sub_df,missing):
    """
    This function keeps only rows where (source_id, member_id, experiment_id)
    exist in df_missing.
    """
    merged = sub_df.merge(missing, on=["source_id", "member_id", "experiment_id"], how="inner")
    return not merged.empty

from typing import Optional, Dict, Callable
from itertools import chain
from pathlib import Path

import re
from datetime import datetime

def extract_time_range_from_url(url):
    # Try 8-digit dates first (YYYYMMDD)
    match = re.search(r"(\d{8})-(\d{8})", url)
    if match:
        start = datetime.strptime(match.group(1), "%Y%m%d")
        end = datetime.strptime(match.group(2), "%Y%m%d")
        return start, end

    # Then try 6-digit dates (YYYYMM)
    match = re.search(r"(\d{6})-(\d{6})", url)
    if match:
        start = datetime.strptime(match.group(1), "%Y%m")
        end = datetime.strptime(match.group(2), "%Y%m")
        return start, end

    # Fallback if no match
    return None, None

def is_time_range_incomplete(url, experiment_id):
    start, end = extract_time_range_from_url(url)
    if not start or not end:
        return True  # Can't determine, assume incomplete
    
    if experiment_id in ["ssp126","ssp245","ssp585"]:
        expected_start = datetime(2015, 1, 16)
        expected_end = datetime(2100, 12, 1)
    elif experiment_id == "historical":
        expected_start = datetime(1850, 1, 16)
        expected_end = datetime(2014, 12, 1)
    else:
        # Add other experiment_id ranges as needed
        return False

    return start > expected_start or end < expected_end

BAD_HOSTS = ["diasjp.net", "esgf-data04.diasjp.net"]

def _norm_url(url):
    url = str(url).replace(".html", "")
    if url.startswith("http://"):
        url = "https://" + url[len("http://"):]
    return url

def _host(url):
    try:
        return urlparse(url).netloc.lower()
    except Exception:
        return ""

def _basename(url):
    try:
        return os.path.basename(urlparse(url).path)
    except Exception:
        return os.path.basename(url)

def _extract_urls_from_file_info(catalog, verbose=True):
    infos = catalog._get_file_info(separator="|", quiet=not verbose)
    urls = []

    for rec in infos:
        if not isinstance(rec, dict):
            continue
        for _, v in rec.items():
            if isinstance(v, str) and ("http://" in v or "https://" in v):
                urls.append(v)
            elif isinstance(v, (list, tuple)):
                for vv in v:
                    if isinstance(vv, str) and ("http://" in vv or "https://" in vv):
                        urls.append(vv)

    # clean + dedupe + drop bad hosts
    out = []
    seen = set()
    for u in urls:
        u = _norm_url(u)
        if u in seen:
            continue
        seen.add(u)
        if any(bad in _host(u) for bad in BAD_HOSTS):
            continue
        out.append(u)

    # prefer fileServer over dodsC
    out = sorted(
        out,
        key=lambda u: (
            0 if "/thredds/fileServer/" in u else 1,
            _host(u),
            u,
        )
    )
    return out

LOCAL_ESGF_CACHE = os.path.expanduser("~/.esgf_manual")

def _local_cache_path(url):
    rel = urlparse(url).path.lstrip("/")
    rel = rel.replace("thredds/fileServer/", "")
    rel = rel.replace("thredds/dodsC/", "")
    return os.path.join(LOCAL_ESGF_CACHE, rel)

def _open_local_dataset(local_path, chunks=None, engine=None, drop_variables=None, **kwargs):
    return xr.open_dataset(
        local_path,
        engine=engine,
        drop_variables=drop_variables,
        chunks=chunks,
        decode_cf=True,
        use_cftime=True,
        **kwargs,
    )

def _download_file_with_retries(url, local_path, n_retries=4, timeout=60, verbose=True):
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    # reuse good cached file
    if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
        if verbose:
            print(f"Using cached file: {local_path}")
        return local_path

    session = requests.Session()

    for attempt in range(1, n_retries + 1):
        tmp_path = local_path + f".part{attempt}"
        try:
            #if verbose:
            #    print(f"Download attempt {attempt}/{n_retries}: {url}")

            with session.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(tmp_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

            if not os.path.exists(tmp_path) or os.path.getsize(tmp_path) == 0:
                raise OSError("Downloaded empty file")

            os.replace(tmp_path, local_path)

            if verbose:
                print(f"Downloaded OK: {local_path}")
            return local_path

        except Exception as e:
            if verbose:
                print(f"Download failed ({attempt}/{n_retries}): {type(e).__name__}: {e}")
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass

            if attempt < n_retries:
                sleep_s = (2 ** (attempt - 1)) + random.uniform(0, 1)
                if verbose:
                    print(f"Sleeping {sleep_s:.1f}s before retry")
                time.sleep(sleep_s)

    raise OSError(f"All download attempts failed for {url}")

def _download_and_open(url, chunks=None, engine=None, drop_variables=None, n_retries=4, verbose=True, **kwargs):
    if "/thredds/fileServer/" in url:
        local_path = _local_cache_path(url)
        local_path = _download_file_with_retries(
            url,
            local_path,
            n_retries=n_retries,
            verbose=verbose,
        )
        return _open_local_dataset(
            local_path,
            chunks=chunks,
            engine=engine,
            drop_variables=drop_variables,
            **kwargs,
        )

    # For dodsC, try a few times too
    last_err = None
    for attempt in range(1, n_retries + 1):
        try:
            if verbose:
                print(f"Remote open attempt {attempt}/{n_retries}: {url}")
            return xr.open_dataset(
                url,
                engine=engine,
                drop_variables=drop_variables,
                chunks=chunks,
                decode_cf=True,
                use_cftime=True,
                **kwargs,
            )
        except Exception as e:
            last_err = e
            if verbose:
                print(f"Remote open failed ({attempt}/{n_retries}): {type(e).__name__}: {e}")
            if attempt < n_retries:
                sleep_s = (2 ** (attempt - 1)) + random.uniform(0, 1)
                time.sleep(sleep_s)

    raise last_err

def load_from_catalog(
    catalog,
    chunks: Optional[Dict[str, int]] = {'time': 200},
    prefer_opendap: bool = False,
    preprocess: Optional[Callable] = None,
    postprocess: Optional[Callable] = None,
    parallel: bool = True,
    engine: Optional[str] = None,
    drop_variables: Optional[list] = None,
    combine: str = 'by_coords',
    concat_dim: Optional[str] = None,
    combine_by_coords_kwargs: Optional[Dict] = None,
    combine_method: str = 'manual',
    max_files: Optional[int] = None,
    verbose: bool = True,
    esgf_url=None,
    n_retries: int = 4,
    **kwargs
):
    try:
        urls = _extract_urls_from_file_info(catalog, verbose=verbose)

        if verbose:
            print(f"Candidate ESGF URLs: {len(urls)}")

        if not urls:
            return None

        if max_files:
            urls = urls[:max_files]

        variant_pat = re.compile(r'_(r\d+i\d+p\d+f\d+)_')
        ym_pat = re.compile(r'_(\d{6})-\d{6}')

        def get_variant(fname):
            m = variant_pat.search(fname)
            return m.group(1) if m else "unknown"

        def get_variable(fname):
            return _basename(fname).split('_')[0]

        def start_ym(fname: str) -> int:
            m = ym_pat.search(fname)
            return int(m.group(1)) if m else 999999

        grouped = defaultdict(list)
        for u in urls:
            grouped[get_variant(u)].append(u)

        datasets = []

        for variant, variant_urls in sorted(grouped.items()):
            by_var = defaultdict(list)
            for u in variant_urls:
                by_var[get_variable(u)].append(u)

            var_dsets = []

            for var, var_urls in sorted(by_var.items()):
                # group replicas by logical file name
                file_groups = defaultdict(list)
                for u in var_urls:
                    file_groups[_basename(u)].append(u)

                time_dsets = []

                for logical_file, candidates in sorted(file_groups.items(), key=lambda kv: start_ym(kv[0])):
                    ds = None
                    for url in candidates:
                        try:
                            #if verbose:
                            #    print("Trying", url)
                            ds = _download_and_open(
                                url,
                                chunks=chunks,
                                engine=engine,
                                drop_variables=drop_variables,
                                n_retries=n_retries,
                                **kwargs,
                            )
                            if preprocess:
                                ds = preprocess(ds)
                            if postprocess:
                                ds = postprocess(ds)
                            break
                        except Exception as e:
                            pass
                            #if verbose:
                            #    print(f"Failed: {type(e).__name__}: {e}")

                    if ds is not None:
                        time_dsets.append(ds)

                if not time_dsets:
                    continue

                var_cat = xr.concat(
                    time_dsets,
                    dim="time",
                    coords="minimal",
                    data_vars="minimal",
                    compat="override",
                    combine_attrs="override",
                )
                var_dsets.append(var_cat)

            if not var_dsets:
                continue

            merged = xr.merge(var_dsets, compat="override")
            datasets.append(merged)

        return datasets if datasets else None

    except Exception as e:
        if verbose:
            print(f"load_from_catalog failed: {type(e).__name__}: {e}")
        return None
        
def load_first_valid_entry(catalog):
    catalog.remove_ensembles()
    urls = _extract_urls_from_file_info(catalog, verbose=False)

    for url in urls:
        try:
            ds = _download_and_open(url)
            return ds
        except Exception:
            pass

    return None

class CMIP6():
    def __init__(self, variable, experiment_id, compare_exp=None, source_id='all', sector_mean=None,sector_sum=None, members=None,time_chunks=None,
                 grid_label=['gn'], new_grid=None, method='conservative_normed',client=False,sic_mask=None,verbose=False,skip_sids=None,table_id=None,
                 esgf_url='https://esgf-node.llnl.gov/esg-search',
                 cat_url="https://cmip6-pds.s3.amazonaws.com/pangeo-cmip6.json",
                 cat_url2="https://storage.googleapis.com/cmip6/cmip6-pgf-ingestion-test/catalog/catalog.json"):
        self.variable = variable
        self.experiment_id = experiment_id
        self.compare_exp = compare_exp
        self.source_id = source_id 
        self.sector_mean = sector_mean
        if self.sector_mean in ['Inner Arctic','IA']:
           self.sector_mean = 'Inner_Arctic'
        self.sector_sum = sector_sum
        if self.sector_sum in ['Inner Arctic','IA']:
           self.sector_sum = 'Inner_Arctic' 
        self.members = members
        self.grid_label = grid_label
        self.new_grid = new_grid
        self.method = method
        self.client = client
        self.sic_mask = sic_mask
        self.verbose = verbose
        self.skip_sids = skip_sids
        if type(self.skip_sids)==str:
            self.skip_sids = [self.skip_sids]
        if type(self.grid_label)==str:
            self.grid_label = [self.grid_label]
        if type(self.experiment_id)==str:
            self.experiment_id = [self.experiment_id]
        if type(self.source_id)==str:
            self.source_id = [self.source_id]
        self.tid = table_id
        #self.load_from_cloud = load_from_cloud
        if self.variable in ['sia','sit','sit_d','sivolume','snt','sifb','sifb_d','sifb_d2','sifb_d3','sic','rhoi','rhoi2'] and self.tid is None:
            self.tid = 'SImon'
        if self.variable in ['tas','ts'] and self.tid is None:
            self.tid = 'Amon'
        if self.variable in ['tos'] and self.tid is None:
            self.tid = 'Omon'
        if time_chunks == None:
            if self.tid in ['SImon','Amon','Omon']:
                self.chunks = 200
            if self.tid=='SIday':
                self.chunks = 50
        if time_chunks != None:
            self.chunks = time_chunks

        self.col = intake.open_esm_datastore(cat_url,registry=dvr)
        self.col2 = intake.open_esm_datastore(cat_url2,registry=dvr)
        self.esgf_url = esgf_url
        self.col_esgf = ESGFCatalog()
            
    def _get_cat(self,vids,tid,grid,detailed=False):
        #we don't really need to reduce source ids at this point unless it is computationally expensive
        #we can search for where members/source_ids match after data is loaded
        if self.compare_exp==None:
            args = {'experiment_id':self.experiment_id,'table_id':tid,'grid_label':[],'source_id':[]}
        else:
            args = {'experiment_id':self.compare_exp,'table_id':tid,'grid_label':[],'source_id':[]}
        if grid!='all':
            args['grid_label'].extend([grid])
        if self.source_id!=['all']:
            args['source_id'].extend(self.source_id)
        args = {k: v for k, v in args.items() if len(v)!=0}
        if self.variable in ['sifb','siitdconc','siitdsnthick','siitdthick','rhoi']:
            c = self.col2
        else:
            c = self.col
        cat_tgt = c.search(**args,variable_id=vids['tgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id)
        cat_awgt = c.search(**args,variable_id=vids['awgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id)
        cat_aswgt = c.search(**args,variable_id=vids['aswgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id)
        #need to add siconc here because sit needs to be converted in 'TaiESM1' by dividing by siconc twice
        #but I just want to add siconc for 'TaiESM1' only (no need to load siconc for all models if sic_mask==None)
        if self.sic_mask==None and self.variable in ['sit'] and (self.source_id[0] in ['all','TaiESM1'] or 'TaiESM1' in self.source_id):
            #args['source_id']=['TaiESM1']
            vids2 = deepcopy(vids)
            [vids2[list(vids2.keys())[x]].append('siconc') for x in [0,1,2] if vids2[list(vids2.keys())[x]][0] is not None]
            cat_tgt2 = c.search(**args,variable_id=vids2['tgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id,source_id='TaiESM1')
            cat_awgt2 = c.search(**args,variable_id=vids2['awgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id,source_id='TaiESM1')
            cat_aswgt2 = c.search(**args,variable_id=vids2['aswgt'],require_all_on=['source_id','member_id']).search(experiment_id=self.experiment_id,source_id='TaiESM1')
            df_tgt = pd.concat([cat_tgt.df,cat_tgt2.df]).drop_duplicates()
            df_awgt = pd.concat([cat_awgt.df,cat_awgt2.df]).drop_duplicates()
            df_aswgt = pd.concat([cat_aswgt.df,cat_aswgt2.df]).drop_duplicates()
            cat_tgt.esmcat._df = df_tgt
            cat_awgt.esmcat._df = df_awgt
            cat_aswgt.esmcat._df = df_aswgt
        cat = {'tgt':cat_tgt,'awgt':cat_awgt,'aswgt':cat_aswgt}

        if self.variable in ['sifb_d','sifb_d2','sifb_d3','sit_d','rhoi','rhoi2']:
            if vids['tgt'][0]!=None:
                vids_tgt_esgf = dvr[self.variable].query['variable_id']
            else:
                vids_tgt_esgf = vids['tgt']
            if vids['awgt'][0]!=None:
                vids_awgt_esgf = dvr[self.variable].query['variable_id']
            else:
                vids_awgt_esgf = vids['awgt']
                
        else: 
            vids_tgt_esgf = vids['tgt']
            vids_awgt_esgf = vids['awgt']

        ##remove_incomplete(complete=require_all) has the same functionality of require_all_on=['source_id','member_id']
        ##it would take too much effort to make this work with compare_exp, so I think I will remove it as an argument–
        ##I don't make use of it anyway!
        try:
            cat_esgf_tgt = ESGFCatalog().search(**args,variable_id=vids_tgt_esgf,quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp=self.experiment_id, var=vids_tgt_esgf))
        except (NoSearchResults, SearchAPIError) as e:
            cat_esgf_tgt = ESGFCatalog().search(
                experiment_id='historical',source_id='CESM2',variable_id='siconc',quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp='ssp245', var=vids_tgt_esgf))
        try:
            cat_esgf_awgt = ESGFCatalog().search(**args,variable_id=vids_awgt_esgf,quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp=self.experiment_id, var=vids_awgt_esgf))
        except (NoSearchResults, SearchAPIError) as e:
            cat_esgf_awgt = ESGFCatalog().search(
                experiment_id='historical',source_id='CESM2',variable_id='siconc',quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp='ssp245', var=vids_tgt_esgf))
        try:
            cat_esgf_aswgt = ESGFCatalog().search(**args,variable_id=vids['aswgt'],quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp=self.experiment_id, var=vids['aswgt']))
        except (NoSearchResults, SearchAPIError) as e:
            cat_esgf_aswgt = ESGFCatalog().search(
                experiment_id='historical',source_id='CESM2',variable_id='siconc',quiet=True).remove_incomplete(
                complete=lambda sub_df: require_all(sub_df, exp='ssp245', var=vids_tgt_esgf))
            
        cat_esgf = {'tgt':cat_esgf_tgt,'awgt':cat_esgf_awgt,'aswgt':cat_esgf_aswgt}

        if vids['awgt'][0]!=None:
            df = pd.concat([cat['tgt'].df,cat['awgt'].df]).drop_duplicates(
                subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True)
            df = df[df['variable_id']==vids['awgt'][0]]
            cat['awgt'].esmcat._df = df
            if vids['aswgt'][0]!=None:
                df1 = pd.concat([cat['tgt'].df,cat['awgt'].df,cat['aswgt'].search(variable_id='siconc').df]).drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True)
                df2 = pd.concat([cat['tgt'].df,cat['awgt'].df,cat['aswgt'].search(variable_id=vids['aswgt'][0]).df]).drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True)
                df1 = df1[df1['variable_id']=='siconc']
                df2 = df2[df2['variable_id']==vids['aswgt'][0]]
                cat['aswgt'].esmcat._df = pd.concat([df1,df2]).reset_index(drop=True)

        cloud_filtered = pd.concat([cat['tgt'].df,cat['awgt'].df,cat['aswgt'].df]).drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True).drop(
            columns=['zstore', 'dcpp_init_year', 'activity_id', 'version','variable_id','institution_id'], errors='ignore')
        esgf_filtered = pd.concat([cat_esgf['tgt'].df,cat_esgf['awgt'].df,cat_esgf['aswgt'].df]).drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True).drop(
            columns=['project','mip_era','activity_drs','id','version','variable_id','institution_id'], errors='ignore')
        esgf_tgt_filtered = cat_esgf['tgt'].df.drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True).drop(
            columns=['project','mip_era','activity_drs','id','version','variable_id','institution_id'], errors='ignore')
        esgf_awgt_filtered = cat_esgf['awgt'].df.drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True).drop(
            columns=['project','mip_era','activity_drs','id','version','variable_id','institution_id'], errors='ignore')
        esgf_aswgt_filtered = cat_esgf['aswgt'].df.drop_duplicates(
                        subset=['source_id','member_id','experiment_id'],keep='first').reset_index(drop=True).drop(
            columns=['project','mip_era','activity_drs','id','version','variable_id','institution_id'], errors='ignore')

        missing_tgt = esgf_filtered.merge(cloud_filtered, 
                                      on=['source_id', 'member_id', 'experiment_id', 'grid_label', 'table_id'], 
                                      how='left', 
                                      indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        missing_awgt = missing_tgt.merge(esgf_tgt_filtered, how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
        missing_aswgt = missing_awgt.merge(esgf_awgt_filtered, how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
                                                                                                                               
        cat_esgf['tgt'].remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=missing_tgt))
        cat_esgf['awgt'].remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=missing_awgt))
        cat_esgf['aswgt'].remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=missing_aswgt))

        return cat,cat_esgf
    
    def _get_vids(self):
        try:
            if self.variable=='sia' and self.sector_sum=='NH':
                vids = {'tgt':['siarean'],'awgt':['siconc'],'aswgt':[None]}
            if self.variable=='sia' and (self.sector_sum=='Inner_Arctic' or isinstance(self.sector_sum, dict)):
                vids = {'tgt':[None],'awgt':['siconc'],'aswgt':[None]}
            if self.variable=='sia' and self.sector_sum=='SH':
                vids = {'tgt':['siareas'],'awgt':['siconc'],'aswgt':[None]}

            if self.variable=='sivolume' and self.sector_sum=='NH':
                vids = {'tgt':['sivoln'],'awgt':['sivol'],'aswgt':['sithick','siconc']}
            if self.variable=='sivolume' and self.sector_sum=='SH':
                vids = {'tgt':['sivols'],'awgt':['sivol'],'aswgt':['sithick','siconc']}
            if self.variable=='sivolume' and (self.sector_sum=='Inner_Arctic' or isinstance(self.sector_sum, dict)):
                vids = {'tgt':[None],'awgt':['sivol'],'aswgt':['sithick','siconc']}
                
            if self.variable=='sic' and (isinstance(self.new_grid,xr.Dataset) or self.new_grid==None):
                vids = {'tgt':['siconc'],'awgt':[None],'aswgt':[None]}
            if self.variable=='sic' and self.sector_sum!=None:
                vids = {'tgt':[None],'awgt':['siconc'],'aswgt':[None]}

            if self.variable=='sit' and (isinstance(self.new_grid,xr.Dataset) or self.new_grid==None):
                vids = {'tgt':['sithick'],'awgt':[None],'aswgt':[None]}
            if self.variable=='sit' and self.sector_mean!=None:
                vids = {'tgt':[None],'awgt':['sithick'],'aswgt':[None]}

            if self.variable=='snt' and (isinstance(self.new_grid,xr.Dataset) or self.new_grid==None):
                vids = {'tgt':['sisnthick'],'awgt':[None],'aswgt':[None]}
            if self.variable=='snt' and self.sector_mean!=None:
                vids = {'tgt':[None],'awgt':['sisnthick'],'aswgt':[None]}

            if self.variable not in ['sia','sivolume','sic','sit','snt'] and (isinstance(self.new_grid,xr.Dataset) or self.new_grid==None):
                vids = {'tgt':[self.variable],'awgt':[None],'aswgt':[None]}
            if self.variable not in ['sia','sivolume','sic','sit','snt'] and self.sector_mean!=None:
                vids = {'tgt':[None],'awgt':[self.variable],'aswgt':[None]}
        
            if self.sic_mask!=None: #and self.variable not in ['sifb_d', 'sifb_d2', 'rhoi', 'rhoi2', 'sit_d']:
                [vids[list(vids.keys())[x]].append('siconc') for x in [0,1,2] if vids[list(vids.keys())[x]][0] is not None]         
            return vids
        except Exception as error:
            if self.verbose==True:
                print("An error occurred:", type(error).__name__, "–", error)
                print('Warning: Variable_id or other attribute error. Check attributes.')
    
    def _convert_to_dataset(self,ds):
        if self.sector_sum != None and not isinstance(self.sector_sum,dict):
            attrs={'sector':self.sector_sum}
        elif self.sector_mean != None and not isinstance(self.sector_mean,dict):
            attrs={'sector':self.sector_mean}
        elif isinstance(self.new_grid,xr.Dataset):
            attrs={'description':'regridded'}
        elif isinstance(self.sector_mean,dict):
            if list(self.sector_mean.keys())[0] == 'Arctic':
                df = NH_seaice_regions
            if list(self.sector_mean.keys())[0] == 'Antarctic':
                df = SH_seaice_regions
            attrs = {'sector':list(self.sector_mean.keys())[0]+[': '+', '.join(m for i,m in enumerate(
                df.Region.where((df.index).isin(list(self.sector_mean.values())[0])).dropna()))][0]}
        elif isinstance(self.sector_sum,dict):
            if list(self.sector_sum.keys())[0] == 'Arctic':
                df = NH_seaice_regions
            if list(self.sector_sum.keys())[0] == 'Antarctic':
                df = SH_seaice_regions
            attrs = {'sector':list(self.sector_sum.keys())[0]+[': '+', '.join(m for i,m in enumerate(
                df.Region.where((df.index).isin(list(self.sector_sum.values())[0])).dropna()))][0]}
        else:
            attrs={'description':'No spatial subsetting or regridding'}
        ds = ds.assign_attrs(attrs)
        ds = ds.to_dataset(name=self.variable)
        if self.variable == 'sifb_d2':
            ds = ds.rename({'sifb_d2':'sifb_d'})
        if self.variable == 'sifb_d3':
            ds = ds.rename({'sifb_d3':'sifb_d'})
        if self.variable == 'sit_d':
            ds = ds.rename({'sit_d':'sit'})
        ds = ds.assign_attrs(attrs)
        return ds
    
    def _fix_grids(self,ds,weight):
        ds=ds.copy()
        if 'lat' in ds.coords and 'lat' in weight.coords and 'y' in ds.coords and 'x' in weight.coords:
            if ds.lat.max().values>90 or 'time' in ds.lat.dims or ds.lat.isnull().sum()>0:
                ds['lat'] = weight.lat
                ds['lon'] = weight.lon
            if weight.y.size != ds.y.size:
                weight = weight.isel(y=slice(0,ds.y.size)
                                 ,x=slice(0,ds.x.size))    
            if ds.isel(y=0).lat.values[0]!=weight.isel(y=0).lat.values[0]:
                if weight.lat.isel(y=0,x=0).values > 0:
                    weight = weight.reindex(y=list(reversed(weight.y)))
                    weight['lat'] = ds.lat
                    weight['lon'] = ds.lon
                    weight['y'] = ds.y
                    weight['x'] = ds.x
                if ds.lat.isel(y=0,x=0).values > 0:
                    ds = ds.reindex(y=list(reversed(weight.y)))
                    ds['lat'] = weight.lat
                    ds['lon'] = weight.lon
                    ds['y'] = weight.y
                    ds['x'] = weight.x
            if ds.y.values[0]!=weight.y.values[0]:
                weight['y'] = ds.y
                weight['x'] = ds.x
            if ds['lat'].mean('x').isnull().any().compute().item():
                ds['lat'] = weight.lat
                ds['lon'] = weight.lon
                ds['y'] = weight.y
                ds['x'] = weight.x
        if self.new_grid!=None:
            if ds.lon_b.max()>90 and 'lat_b' in weight.coords:
                ds['lat_b'] = weight.lat_b
                ds['lon_b'] = weight.lon_b
                ds['y_b'] = weight.y_b
                ds['x_b'] = weight.x_b
                ds['lat'] = weight.lat
                ds['lon'] = weight.lon
                ds['y'] = weight.y
                ds['x'] = weight.x
        return ds,weight

    def _remove_vars(self, ds):
        ds = ds.copy()
        if 'dcpp_init_year' in ds.dims:
            ds = ds.isel(dcpp_init_year=0)  
        if 'nodes' in ds.dims:
            ds=ds.isel(nodes=0)
        if self.sector_mean is not None or self.sector_sum is not None:
            if 'x' in ds.dims and ds.x.size==1:
                ds=ds.isel(x=0,y=0).drop_vars(['x','y','lat','lon'], errors="ignore")
            if 'vertex' in ds.dims:
                ds=ds.isel(vertex=0)
            if 'bnds' in ds.dims:
                ds=ds.isel(bnds=0)
        ds = ds.drop_vars(['sector','dcpp_init_year','nodes','type','time_bounds','vertex','bnds','iceband_bnds'], errors="ignore") 
        #else:
        #    ds = ds.drop_vars(['sector','dcpp_init_year','nodes','type','time_bounds'], errors="ignore") 
        return ds

    def _postprocessing(self, ds):
        ds = ds.copy()
        if 'time' in ds.variables:
            ds = convert_time(ds)
            #ds = convert_time2(ds)
            ds = update_member_id(ds)
            #if 'time' in ds.dims and all(var.chunks is None for var in ds.data_vars.values()):
            #ds = set_chunks(ds,self.chunks)
            ds = ds.drop_vars(['area','sector'],errors='ignore')
            #if 'units' in ds.data_vars.dtypes:
            try:
                if ds[next(iter(ds.data_vars.dtypes))].units =='1e3 km3':
                    if ds[next(iter(ds.data_vars.dtypes))].max() > 1000:
                        ds = ds/1000.
                    if ds[next(iter(ds.data_vars.dtypes))].max() > 1e8:
                        ds = ds/1e9
                    if 'y' in ds.dims:
                        ds = ds.max(['x','y'])
            except Exception as error:         
                pass
        return ds 
        
    def _get_area(self,sid,grid):
        #no areacello for 'UKESM1-1-LL', so use 'UKESM1-0-LL' instead (same grid it appears)
        if sid == 'UKESM1-1-LL':
            sid = 'UKESM1-0-LL'
        #the only areacello for NESM3 has wrong resolution for some reason
        #CESM2 appears to have the same grid, with maybe slightly different ocean masking
        #f sid == 'NESM3':
        #    sid = 'CESM2'
        area_subset = self.col.search(experiment_id=self.experiment_id, variable_id=['areacello'],
                                 source_id=[sid],table_id=['Ofx'], grid_label=grid)
        if area_subset.df.source_id.size == 0:
            area_subset = self.col.search(variable_id=['areacello'],
                                     source_id=[sid],table_id=['Ofx'], grid_label=grid) 
        if self.variable in ['siflswutop','siflswdtop','siflswdbot'] and sid in ['ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-LR']:
            area_subset = self.col.search(variable_id=['areacella'],
                                source_id=[sid],table_id=['fx'], grid_label=grid)
        df = area_subset.df.groupby(['source_id']).first().reset_index()
        area_subset.esmcat._df = df
        area_dict = area_subset.to_dataset_dict(aggregate=True,xarray_open_kwargs={"consolidated": True, 'decode_times':False}
                                                ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False)
        if bool(area_dict):
            area = area_dict[list(area_dict.keys())[0]].squeeze()
            area = area.where(area<1e35)
            area = area.drop_vars(['member_id','dcpp_init_year'], errors="ignore")
            if 'areacella' in area.data_vars:
                area = area.rename({'areacella':'areacello'})
            return area

        try:
            area_subset = ESGFCatalog().search(source_id=sid,variable_id='areacello',grid_label=grid,quiet=True)
            area = load_first_valid_entry(area_subset)
            if bool(area):
                area = complete_preprocessing(area)
                area = area.where(area<1e35)
                area = area.drop_vars(['member_id','dcpp_init_year'], errors="ignore")
                return area
        except NoSearchResults:
            pass  # skip silently if no data
    
    def _get_ocean_mask(self,sid,grid):
#        if sid == 'UKESM1-1-LL':
#            sid = 'UKESM1-0-LL'
#        if sid == 'NESM3':
#            sid = 'CESM2'
        ocean_subset = self.col.search(experiment_id=self.experiment_id, variable_id=['sftof'],
                                  source_id=[sid],table_id=['Ofx'], grid_label=grid)
        if ocean_subset.df.source_id.size == 0:
            ocean_subset = self.col.search(variable_id=['sftof'],
                                      source_id=[sid],table_id=['Ofx'], grid_label=grid)
        if self.variable in ['siflswutop','siflswdtop','siflswdbot'] and sid in ['ACCESS-CM2','MRI-ESM2-0','MPI-ESM1-2-LR']:
            ocean_subset = self.col.search(variable_id=['sftlf'],
                                source_id=[sid],table_id=['fx'], grid_label=grid)
        df = ocean_subset.df.groupby(['source_id']).first().reset_index()
        ocean_subset.esmcat._df = df
        ocean_dict = ocean_subset.to_dataset_dict(aggregate=True,xarray_open_kwargs={"consolidated": True, 'decode_times':False}
                                                  ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False)
        if bool(ocean_dict):
            ocean = ocean_dict[list(ocean_dict.keys())[0]].squeeze()
            if 'sftlf' in ocean.data_vars:
                ocean.rename({'sftlf':'sftof'})
                ocean = (100-ocean)*100
            return ocean.drop_vars(['member_id','dcpp_init_year'], errors="ignore")

        try:
            ocean_subset = ESGFCatalog().search(source_id=sid,variable_id='sftof',grid_label=grid,quiet=True)
            ocean = load_first_valid_entry(ocean_subset)
            if bool(ocean):
                ocean = complete_preprocessing(ocean)
                return ocean.drop_vars(['member_id','dcpp_init_year'], errors="ignore")
        except NoSearchResults:
            pass  # skip silently if no data
    
    def _spatial_average(self,ds,weight=None,sic=None):
        ds = ds.copy()
        if 'lat' in ds.coords and 'y' in ds.dims:
            if self.sector_mean == 'NH':
                ds_subset = ds.where(ds.lat>0)
                lat = ds.lat.where(ds.lat>0)
            elif self.sector_mean == 'SH':
                ds_subset = ds.where(ds.lat<0)
                lat = ds.lat.where(ds.lat<0)
            elif self.sector_mean == 'Arctic':
                ds_subset = ds.where(ds.lat>60)
                lat = ds.lat.where(ds.lat>60)
            elif self.sector_mean == 'Antarctic':
                ds_subset = ds.where(ds.lat<-60)
                lat = ds.lat.where(ds.lat<-60)
            elif self.sector_mean == 'Inner_Arctic':
                df = NH_seaice_regions
                mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
                ds_subset = ds.where(mask.isin([0,1,2,3,4,5,6]))
                lat = ds.lat.where(mask.isin([0,1,2,3,4,5,6]))
            elif isinstance(self.sector_mean, dict): 
                if list(self.sector_mean.keys())[0] == 'Arctic':
                    df = NH_seaice_regions
                    mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
                    ds_subset = ds.where(mask.isin(list(self.sector_mean.values())[0]))
                    lat = ds.lat.where(mask.isin(list(self.sector_mean.values())[0]))
            elif isinstance(self.sector_mean, dict): 
                if list(self.sector_mean.keys())[0] == 'Antarctic':
                    df = SH_seaice_regions
                    mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
                    ds_subset = ds.where(mask.isin(list(self.sector_mean.values())[0]))
                    lat = ds.lat.where(mask.isin(list(self.sector_mean.values())[0]))
            else: 
                ds_subset = ds
            if self.sic_mask!=None:
                ds_subset = ds_subset.where(sic>self.sic_mask)
            if not isinstance(weight,xr.DataArray):
                #lat is already a 2D field, so each point is weighted without having to broadcast!
                #if data is on a recticular grid, then we can just weight is by the cosine of lat
                lat = lat.where(ds_subset.notnull())
                weight=np.cos(np.deg2rad(lat))/np.cos(np.deg2rad(lat)).mean('y')
                ds_mean = (ds_subset*weight).mean(['x','y'])
            else:
                #if data is on a curvilinear grid, we need to weight it by the grid-cell area
                #here, we just compute the average over the area where there is data or where 
                #sea ice variables are not 0
                ds_subset,weight = self._fix_grids(ds_subset,weight)
                #should we add sifb here, or are we allowing negative freeboard values?
                if self.variable in ['sit','snt','sic']: 
                    ds_mean = ((ds_subset*weight).sum(['x', 'y'],min_count=1)
                               /weight.where(np.logical_and(ds_subset.notnull(),ds_subset>0)).sum(['x','y']))
                else:
                    ds_mean = ((ds_subset*weight).sum(['x', 'y'],min_count=1)
                               /weight.where(ds_subset.notnull()).sum(['x','y']))
        return ds_mean
            
    def _weighted_sum(self,ds,weight,vid,calc=True):
        ds=ds.copy()
        ds,weight = self._fix_grids(ds,weight)
        if self.sector_sum == 'NH':
            ds_subset = ds.where(ds.lat>0)
        if self.sector_sum == 'SH':
            ds_subset = ds.where(ds.lat<0)
        if self.sector_sum == 'Inner_Arctic':
            df = NH_seaice_regions
            mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
            ds_subset = ds.where(mask.isin([0,1,2,3,4,5,6]))
        elif isinstance(self.sector_sum, dict): 
            if list(self.sector_sum.keys())[0] == 'Arctic':
                df = NH_seaice_regions
                mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
                ds_subset = ds.where(mask.isin(list(self.sector_sum.values())[0]))
            if  list(self.sector_sum.keys())[0] == 'Antarctic':  
                df = SH_seaice_regions
                mask = regionmask.mask_geopandas(df, ds.lon, ds.lat,overlap=False)
                ds_subset = ds.where(mask.isin(list(self.sector_sum.values())[0]))
        if vid == 'siconc':
            ds_subset = ds_subset/100
        if calc==False:
            #weight = (ds_subset*weight.fillna(0))
            weight = (ds_subset*weight)
            return weight
        if calc==True:
            #weighted_sum = (ds_subset*weight.fillna(0)).sum(['x', 'y'],min_count=1)*10**-12
            weighted_sum = (ds_subset*weight).sum(['x', 'y'],min_count=1)*10**-12
            return weighted_sum
    
    def _regrid_old(self,ds,sid,new_grid=None):
        unstr_sids = ['ICON-ESM-LR','AWI-CM-1-1-MR', 'AWI-ESM-1-1-LR']
        if sid not in unstr_sids:
            if ds.lon.isel(x=0).mean()==360 and ds.lon.isel(x=1).mean()==1.:
                ds = ds.roll(x=-1,roll_coords=True)
        if self.new_grid == None:
            ds_out=new_grid
        else:
            ds_out=self.new_grid
        #if 'member_id' in ds.mask.dims:
        #    ds['mask'] = ds.mask.isel(member_id=0)
        if self.method=='bilinear' and sid not in unstr_sids:
            ds['mask'] = ds.mask.where(ds.mask>0)
            ds['mask'].loc[dict(y=ds['mask'].y[-1])] = ds['mask'].loc[dict(y=ds['mask'].y[-1])].fillna(0)
            extrap_method='inverse_dist'
        #   extrap_method="nearest_s2d"
            locstream_in=False
            skipna=False
            method = self.method
        elif self.method in ['conservative_normed','conservative'] and sid not in unstr_sids:
            if (ds.lon_b==180).sum()>ds.y.size:
                #for some reason land masks as 180 here so we need to
                #redo lon_b and lat_b with unmasked lon/lat coords for accurate regridding
                ds = add_corners(ds.drop_vars(['lon_b','lat_b']))
            locstream_in=False
            skipna=True
            extrap_method=None
            method = self.method
        elif sid in unstr_sids:
            #unstructered grids can only be regridded using nearest_s2d
            method = 'nearest_s2d'
            locstream_in=True
            skipna=True
            extrap_method=None
        else:
            locstream_in=False
            extrap_method=None
            skipna=True
            method =self.method
        if ds.mask.shape != ds.lon.shape:
            ds['mask'] = ds.mask.transpose('x','y')
        regridder = xe.Regridder(ds, ds_out, method
                                 ,ignore_degenerate=True,extrap_method=extrap_method,periodic=True,locstream_in=locstream_in)
        ds_new = regridder(ds,skipna=skipna)
        #return ds_new.where(ds_new.mask).chunk({'y':ds_new.y.size,'x':ds_new.x.size})
        return ds_new.where(ds_new.mask)

    def _regrid(self, ds, sid, new_grid=None):
        unstr_sids = ['ICON-ESM-LR','AWI-CM-1-1-MR', 'AWI-ESM-1-1-LR']
        is_tas_like = self.variable in ["tas", "ts"]
    
        # Fix dateline wrap for structured grids
        if sid not in unstr_sids:
            if ds.lon.isel(x=0).mean() == 360 and ds.lon.isel(x=1).mean() == 1.:
                ds = ds.roll(x=-1, roll_coords=True)        
    
        ds_out = new_grid if self.new_grid is None else self.new_grid
    
        # -------------------------
        # Choose method parameters
        # -------------------------
        if is_tas_like:
            # tas/ts: bilinear is the right default; no mask required
            method = "bilinear"
            locstream_in = False
            skipna = True
            extrap_method = "nearest_s2d"  # or None
    
        elif self.method == 'bilinear' and sid not in unstr_sids:
            # sea-ice/ocean vars with a mask
            if "mask" in ds:
                ds["mask"] = ds["mask"].where(ds["mask"] > 0)
                ds["mask"].loc[dict(y=ds["mask"].y[-1])] = ds["mask"].loc[dict(y=ds["mask"].y[-1])].fillna(0)
    
            method = self.method
            locstream_in = False
            skipna = True
            extrap_method = "nearest_s2d"  # was inverse_dist; nearest is safer
    
        elif self.method in ['conservative_normed','conservative'] and sid not in unstr_sids:
            #if ("lon_b" in ds.coords) and ((ds.lon_b == 180).sum() > ds.y.size):
            #    ds = add_corners(ds.drop_vars(['lon_b','lat_b'], errors="ignore"))
            # More robust seam detection
            if ("lon_b" in ds.coords) and ((np.abs(ds.lon_b - 180) < 0.5).sum() > ds.y.size):
                # Find seam location robustly and roll to edge before add_corners
                lon_diff = np.abs(ds.lon.diff('x'))
                seam_x = int(lon_diff.mean('y').argmax('x').values)
                ds = ds.roll(x=-(seam_x + 1), roll_coords=True)
                ds = add_corners(ds.drop_vars(['lon_b','lat_b'], errors="ignore"))
    
            method = self.method
            locstream_in = False
            skipna = True
            extrap_method = None
    
        elif sid in unstr_sids:
            method = 'nearest_s2d'
            locstream_in = True
            skipna = True
            extrap_method = None
    
        else:
            method = self.method
            locstream_in = False
            skipna = True
            extrap_method = None
    
        # -------------------------
        # Only handle mask if present and needed
        # -------------------------
        if (not is_tas_like) and ("mask" in ds):
            if ds["mask"].shape != ds.lon.shape:
                ds["mask"] = ds["mask"].transpose("x", "y")
    
        regridder = xe.Regridder(
            ds, ds_out, method,
            ignore_degenerate=True,
            extrap_method=extrap_method,
            periodic=True,
            locstream_in=locstream_in
        )
    
        ds_new = regridder(ds, skipna=skipna)
    
        # Apply mask only for sea-ice/ocean style vars
        if (not is_tas_like) and ("mask" in ds_new):
            return ds_new.where(ds_new["mask"])
        return ds_new    

    def _load_data(self,vids,cat,cat_esgf,tid,grid):
        if self.variable in ['sifb_d','sifb_d2','sifb_d3','sit_d','rhoi','rhoi2']:
            agg=True
            xr_dict = {'coords':'minimal','data_vars':'minimal','compat':'override','combine_attrs':'override'}
        else:
            agg=False
            xr_dict = {}
        datasets = []
        #I used to skip these, but I keep them despite some regridding issues
        #skip_sids = ['BCC-CSM2-MR','BCC-ESM1','CAS-ESM2-0']
        
        #['ICON-ESM-LR','AWI-CM-1-1-MR', 'AWI-ESM-1-1-LR'] all are on an unstructured grid and need to be regridded in a specific way
        unstr_sids = ['ICON-ESM-LR','AWI-CM-1-1-MR', 'AWI-ESM-1-1-LR']
        if vids['tgt'][0] is not None:
            #if self.load_from_cloud == True:
            source_ids_cloud = cat['tgt'].df['source_id'].unique()
            source_ids_esgf = cat_esgf['tgt'].df['source_id'].unique()
            source_ids = np.union1d(source_ids_cloud,source_ids_esgf)
            if self.verbose==True:
                print('{} loaded from the cloud: {}'.format(vids['tgt'],source_ids_cloud.tolist()))
                print('{} loaded via OpenDap: {}'.format(vids['tgt'],source_ids_esgf.tolist()))
                
            for sid in (pbar := tqdm(source_ids,leave=False)):
                pbar.set_description(f"{self.variable}: {sid}")
            #for sid in source_ids:
                #if sid in skip_sids:
                #    continue
                def add_tgt(ds):
                    if sid=='TaiESM1' and self.variable in ['sit']:
                        ds['sithick']=ds['sithick']/ds['siconc']/(ds['siconc']/100)                  
                    if self.variable in ['sit','snt','sifb','sifb_d','sifb_d2','sifb_d3','sit_d']:
                        ds[vids['tgt'][0]] = ds[vids['tgt'][0]].where(lambda x:np.abs(x)<10)
                    #ds = ds.drop_vars(['area','sector'],errors='ignore')
                    if self.sic_mask!=None and self.new_grid==None:
                        sic = ds['siconc']/100
                        ds = ds.where(sic>self.sic_mask)
                    if self.sector_mean is not None:
                        if 'lat' in ds.coords and 'y' in ds.dims:
                            ds_mean = self._spatial_average(ds[vids['tgt'][0]])
                            ds_mean = self._remove_vars(ds_mean)
                            #if sid=='TaiESM1' and self.variable == 'sifb_d':
                            #    ds_mean['sifb_d'] = ds_mean.sifb_d/100
                            datasets.append(ds_mean)
                        else:
                            pass
                    elif self.new_grid is not None:
                        if self.variable not in ['tas','ts']:
                            try:
                                if 'lon_bounds' in ds.coords or 'lon_b' not in ds.coords:
                                    ds = add_corners(ds.drop_vars(['lon_bounds','lat_bounds','vertex','bounds'],errors="ignore"))
                                area = self._get_area(sid,grid)
                                #if area is not None:
                                #    ds,area = self._fix_grids(ds,area)
                                try:
                                    area.compute()
                                except Exception:
                                    area = calc_areacello(ds)
                                if area is None or sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1','NESM3']:
                                    area = calc_areacello(ds)
                                ds,area = self._fix_grids(ds,area)
                                ds = ds.assign_coords({'areacello':area.areacello})
                                ocean = self._get_ocean_mask(sid,grid)
                                ds,ocean = self._fix_grids(ds,ocean)
                                ds['mask']=(ocean.sftof/100)
                                #ds['x'] = np.arange(0,len(ds.x))
                                #ds['y'] = np.arange(0,len(ds.y))
                            except (IndexError,AttributeError,ValueError):
                                if self.sic_mask!=None:
                                    sic = ds['siconc']/100
                                    ds['mask'] = sic.isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                    ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                                else:
                                    if 'iceband' in ds.coords:
                                        ds['mask'] = ds[vids['tgt'][0]].isel(time=slice(0,120),iceband=0).mean(['time']).isel(member_id=0).notnull().squeeze(
                                        ).drop_vars(['experiment_id','time','member_id','dcpp_init_year','iceband'],errors="ignore")
                                    else:
                                        ds['mask'] = ds[vids['tgt'][0]].isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                        ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                        try:
                            ds_new = self._regrid(ds,sid)
                            if self.sic_mask!=None and isinstance(self.new_grid,xr.Dataset):
                                sic = ds_new['siconc']/100
                                #if sid=='UKESM1-0-LL':
                                #    ds_new = interpolate_na(ds_new, ["y", "x"], method="linear")
                                ds_new = ds_new.where(sic>self.sic_mask)
                            ds_new = self._remove_vars(ds_new)
                            datasets.append(ds_new[vids['tgt'][0]])
                        except Exception as error:
                            if self.verbose==True:
                                print("An error occurred:", type(error).__name__, "–", error)
                                print('Warning: Regridding failed for {} (grid={})'.format(sid,grid))
                    elif self.variable in ['sia','sivolume']:
                        ds = self._remove_vars(ds)
                        datasets.append(ds[vids['tgt'][0]])
                    else:
                        if self.variable not in ['tas','ts'] and sid not in unstr_sids:
                            try:
                                area = self._get_area(sid,grid)
                                #if area is not None:
                                #    ds,area = self._fix_grids(ds,area)
                                try:
                                    area.compute()
                                except Exception:
                                    area = calc_areacello(ds)
                                if area is None or sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1','NESM3']:
                                    area = calc_areacello(ds)
                                ds,area = self._fix_grids(ds,area)
                                ds = ds.assign_coords({'areacello':area.areacello})
                            except (IndexError,AttributeError,ValueError):
                                if self.verbose==True:
                                    print('Warning: areacello from {} did not load (grid={})'.format(sid,grid))
                                ds = ds.assign_attrs(areacello='not found for this model')
                            try:
                                ocean = self._get_ocean_mask(sid,grid)
                                mask = ocean.drop_vars(['member_id','dcpp_init_year'], errors="ignore").sftof
                            except (IndexError,AttributeError,ValueError):
                                if self.sic_mask!=None:
                                    mask = sic.isel(time=0,member_id=0).notnull().squeeze(
                                    ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                                else:
                                    if 'iceband' in ds.coords:
                                        mask = ds[vids['tgt'][0]].isel(time=slice(0,120),iceband=0).mean(['time']).isel(member_id=0).notnull().squeeze(
                                        ).drop_vars(['experiment_id','time','member_id','dcpp_init_year','iceband'],errors="ignore")
                                    else:
                                        mask = ds[vids['tgt'][0]].isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                        ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                            #if sid=='ACCESS-CM2' and self.variable in ['sifb']:
                            #    ds['sifb'] = ds['sifb']/(ds['siconc']/100)
                            ds = ds.assign_coords({'mask':mask})
                        if sid in unstr_sids:
                            try:
                                ocean = self._get_ocean_mask(sid,grid)
                                ds['mask']=(ocean.sftof/100)
                            except (IndexError,AttributeError,ValueError):
                                if 'iceband' in ds.coords:
                                    ds['mask'] = ds[vids['tgt'][0]].isel(time=slice(0,120),iceband=0).mean(['time']).isel(member_id=0).notnull().squeeze(
                                    ).drop_vars(['experiment_id','time','member_id','dcpp_init_year','iceband'],errors="ignore")
                                else:
                                    ds['mask'] = ds[vids['tgt'][0]].isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                        ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore") 
                            ds = self._regrid(ds,sid,grid_CESM2)
                        datasets.append(ds[vids['tgt'][0]])
    
                if sid in source_ids_cloud:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat['tgt'].search(source_id=sid)
                        if self.members=='first':
                            if subset.search(member_id='r1i1p1f1').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f1')
                            elif subset.search(member_id='r1i1p1f2').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f2')
                            elif subset.search(member_id='r1i1p1f3').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f3')
                            else:
                                subset = subset.search(member_id='r1i1p2f1')
                            agg=True
                        try:
                            ds_dict = subset.to_dataset_dict(aggregate=agg
                                                             ,xarray_open_kwargs={"consolidated": True,'decode_times': True,'use_cftime':True, 'chunks': {'time': self.chunks}}
                                                                ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False
                                                        ,xarray_combine_by_coords_kwargs=xr_dict)
                        except Exception as error:
                            ds_dict = subset.to_dataset_dict(aggregate=agg
                                                             ,xarray_open_kwargs={"consolidated": True,'decode_times': True,'use_cftime':False,'chunks': {'time': self.chunks}}
                                                                ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False
                                                        ,xarray_combine_by_coords_kwargs=xr_dict)
                        #print(ds_dict)
                        dsets = list(map(self._postprocessing,ds_dict.values()))
                        #print(dsets)
                        dsets = list(map(self._remove_vars,dsets))
                        #print(dsets)
                        #dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                        #        ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #        for exp in self.experiment_id]
                        #dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                        #        ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #        for exp in self.experiment_id]
                        dsets = [xr.concat([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                                ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override',dim='member_id',join='outer') 
                                for exp in self.experiment_id]
                        #print(dsets)
                        [add_tgt(ds) for ds in dsets if len(ds.data_vars)!=0]
                    
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)
                            print('Warning:{} from {} did not load (grid={})'.format(vids['tgt'],sid,grid))
                        continue
                        
                if sid in source_ids_esgf:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat_esgf['tgt'].clone().search(**ChainMap({'source_id':sid}, cat_esgf['tgt'].last_search),quiet=True
                                                               ).remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=cat_esgf['tgt'].df))
                        if self.members=='first':
                            subset.remove_ensembles()
                    
                        #ds_dict = subset.to_dataset_dict(prefer_streaming=True, add_measures=False, quiet=True)
                        #dsets = list(map(complete_preprocessing,ds_dict.values()))
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = load_from_catalog(
                            catalog=subset,
                            chunks = {'time':self.chunks},
                            preprocess=complete_preprocessing,
                            postprocess=self._postprocessing,
                            prefer_opendap=False,
                            combine_method='manual',
                            esgf_url = self.esgf_url
                        )
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = list(map(self._remove_vars,dsets))
                        #try:
                        #    dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                        #            ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #            for exp in self.experiment_id]
                        #except Exception as error:
                        #    dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                        #        ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #        for exp in self.experiment_id]
                        dsets = [xr.concat([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                                    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                                    for exp in self.experiment_id]
                        if self.variable in ['sifb_d','sifb_d2','sifb_d3','sit_d','rhoi','rhoi2']:
                        #    if any(var not in dsets[0].data_vars for var in vars): ###no longer needed (I think), but needs to be revised if needed
                        #        if self.verbose==True:
                        #            print('Warning: {} not loaded from {}'.format(vids['tgt'],sid))
                        #        continue
                            dsets = list(map(dvr[self.variable].func,dsets))                    
                        [add_tgt(ds) for ds in dsets if len(ds.data_vars)!=0]
                    
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)
                            print('Warning:{} from {} did not load (grid={})'.format(vids['tgt'],sid,grid))
                        continue
                
        if vids['awgt'][0] is not None:
            source_ids_cloud = cat['awgt'].df['source_id'].unique()
            source_ids_esgf = cat_esgf['awgt'].df['source_id'].unique()
            source_ids = np.union1d(source_ids_cloud,source_ids_esgf)
            if self.verbose==True:
                print('{} loaded from the cloud: {}'.format(vids['awgt'],source_ids_cloud.tolist()))
                print('{} loaded via OpenDap: {}'.format(vids['awgt'],source_ids_esgf.tolist()))  
            
            for sid in (pbar := tqdm(source_ids,leave=False)):
                pbar.set_description(f"{self.variable}: {sid}")
            #for sid in source_ids:
                #if sid in skip_sids:
                #    continue

                def add_awgt(ds):
                    if sid=='TaiESM1' and self.variable in ['sit']:
                        ds['sithick']=ds['sithick']/ds['siconc']/(ds['siconc']/100)
                    if self.variable in ['sit','snt','sifb','sifb_d','sifb_d2','sifb_d3','sit_d']:
                        ds[vids['awgt'][0]] = ds[vids['awgt'][0]].where(lambda x:np.abs(x)<10)
                    if sid in unstr_sids:
                        try:
                            ocean = self._get_ocean_mask(sid,grid)
                            ds,ocean = self._fix_grids(ds,ocean)
                            ds['mask']=(ocean.sftof/100)
                        except (IndexError,AttributeError,ValueError):
                            #ds['mask'] = ds[vids['awgt'][0]].mean(['time']).isel(member_id=0).notnull().squeeze(
                            #).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                            if 'iceband' in ds.coords:
                                ds['mask'] = ds[vids['awgt'][0]].isel(time=slice(0,120),iceband=0).mean(['time']).isel(member_id=0).notnull().squeeze(
                                ).drop_vars(['experiment_id','time','member_id','dcpp_init_year','iceband'],errors="ignore")
                            else:
                                ds['mask'] = ds[vids['awgt'][0]].isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore") 
                        ds = self._regrid(ds,sid,new_grid=new_grid)
                    if self.sector_mean is not None:
                        if self.sic_mask==None:
                            sic=None
                        else:
                            sic = ds['siconc']/100
                        ds_mean = self._spatial_average(ds=ds[vids['awgt'][0]],weight=area,sic=sic)
                        ds_mean = self._remove_vars(ds_mean)
                        datasets.append(ds_mean)
                    if self.sector_sum is not None:   
                        ds_sum = self._weighted_sum(ds=ds[vids['awgt'][0]],weight=area,vid=vids['awgt'][0])
                        ds_sum = self._remove_vars(ds_sum)
                        datasets.append(ds_sum) 
                    
                if sid in source_ids_cloud:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat['awgt'].search(source_id=sid)
                        if self.members=='first':
                            if subset.search(member_id='r1i1p1f1').df.member_id.size==1:
                                subset = subset.search(member_id='r1i1p1f1')
                            elif subset.search(member_id='r1i1p1f2').df.member_id.size==1:
                                subset = subset.search(member_id='r1i1p1f2')
                            elif subset.search(member_id='r1i1p1f3').df.member_id.size==1:
                                subset = subset.search(member_id='r1i1p1f3')
                            else:
                                subset = subset.search(member_id='r1i1p2f1')
                            agg=True
                        try:
                            ds_dict = subset.to_dataset_dict(aggregate=agg
                                                             ,xarray_open_kwargs={"consolidated": True,'decode_times': True,'use_cftime':True,'chunks': {'time': self.chunks}}
                                                                    ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False
                                                        ,xarray_combine_by_coords_kwargs=xr_dict)
                        except Exception as error:
                            ds_dict = subset.to_dataset_dict(aggregate=agg
                                                             ,xarray_open_kwargs={"consolidated": True,'decode_times': True,'use_cftime':False,'chunks': {'time': self.chunks}}
                                                                ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False
                                                        ,xarray_combine_by_coords_kwargs=xr_dict)
                        dsets = list(map(self._postprocessing,ds_dict.values()))
                        dsets = list(map(self._remove_vars,dsets))
                        #dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                        #    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #    for exp in self.experiment_id]
                        #dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                        #    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #    for exp in self.experiment_id]
                        dsets = [xr.concat([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                            ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override',dim='member_id',join='outer') 
                            for exp in self.experiment_id]                        
                        [add_awgt(ds) for ds in dsets if len(ds.data_vars)!=0 and ('x' and 'y' in ds.coords or sid in unstr_sids)]
                        
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)
                            print('Warning: {} from {} did not load (grid={})'.format(vids['awgt'],sid,grid))
                        continue

                if sid in source_ids_esgf:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat_esgf['awgt'].clone().search(**ChainMap({'source_id':sid}, cat_esgf['awgt'].last_search),quiet=True
                                                               ).remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=cat_esgf['awgt'].df))
                        if self.members=='first':
                            subset.remove_ensembles()
                    
                        #ds_dict = subset.to_dataset_dict(prefer_streaming=True, add_measures=False, quiet=True)
                        #dsets = list(map(complete_preprocessing,ds_dict.values()))
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = load_from_catalog(
                            catalog=subset,
                            chunks = {'time':self.chunks},
                            preprocess=complete_preprocessing,
                            postprocess=self._postprocessing,
                            prefer_opendap=False,
                            combine_method='manual')
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = list(map(self._remove_vars,dsets))
                        try:
                            dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                                    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                                    for exp in self.experiment_id]
                        except Exception as error:
                            dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                                ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                                for exp in self.experiment_id]
                        if self.variable in ['sifb_d','sifb_d2','sifb_d3','sithick_d','rhoi','rhoi2']:
                        #    if any(var not in dsets[0].data_vars for var in vars): ###no longer needed (I think), but needs to be revised if needed
                        #        if self.verbose==True:
                        #            print('Warning: {} not loaded from {}'.format(vids['awgt'],sid))
                        #        continue
                            dsets = list(map(dvr[self.variable].func,dsets))
                        [add_awgt(ds) for ds in dsets if len(ds.data_vars)!=0 and ('x' and 'y' in ds.coords or sid in unstr_sids)]
                        
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)
                            print('Warning:{} from {} did not load (grid={})'.format(vids['awgt'],sid,grid))
                        continue    
                try:
                    if sid in unstr_sids:
                        new_grid = grid_CESM2
                        area = grid_CESM2.areacello
                    else:
                        area = self._get_area(sid,grid)
                        try:
                            area.compute()
                        except Exception:
                            area = calc_areacello(dsets[0])
                        if area is None or sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1']:
                            #if sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1']:
                            #    dset,area = self._fix_grids(dsets[0],area.areacello)
                            area = calc_areacello(dsets[0])
                        area = area.areacello
                except (IndexError,AttributeError,ValueError):
                    if self.verbose==True:
                        print('Warning: area from {} not found (grid={})'.format(sid,grid))
                    continue
                        
        if vids['aswgt'][0] is not None:
            source_ids_cloud = cat['aswgt'].df['source_id'].unique()
            source_ids_esgf = cat_esgf['aswgt'].df['source_id'].unique()
            source_ids = np.union1d(source_ids_cloud,source_ids_esgf)
            if self.verbose==True:
                print('{} loaded from the cloud: {}'.format(vids['aswgt'],source_ids_cloud.tolist()))
                print('{} loaded via OpenDap: {}'.format(vids['aswgt'],source_ids_esgf.tolist()))  
                
            for sid in (pbar := tqdm(source_ids,leave=False)):
                pbar.set_description(f"{self.variable}: {sid}")
                #if sid in skip_sids:
                #    continue

                def add_aswgt(ds):
                    if sid=='TaiESM1':
                        ds['sithick']=ds['sithick']/ds['siconc']/(ds['siconc']/100)
                    if sid in unstr_sids:
                        try:
                            ocean = self._get_ocean_mask(sid,grid)
                            ds,ocean = self._fix_grids(ds,ocean)
                            ds['mask']=(ocean.sftof/100)
                        except (IndexError,AttributeError,ValueError):
                            #ds['mask'] = ds[vids['aswgt'][0]].mean(['time']).isel(member_id=0).notnull().squeeze(
                            #).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")
                            if 'iceband' in ds.coords:
                                ds['mask'] = ds[vids['aswgt'][0]].isel(time=slice(0,120),iceband=0).mean(['time']).isel(member_id=0).notnull().squeeze(
                                ).drop_vars(['experiment_id','time','member_id','dcpp_init_year','iceband'],errors="ignore")
                            else:
                                ds['mask'] = ds[vids['aswgt'][0]].isel(time=slice(0,120)).mean(['time']).isel(member_id=0).notnull().squeeze(
                                ).drop_vars(['experiment_id','time','member_id','dcpp_init_year'],errors="ignore")

                        ds = self._regrid(ds,sid,new_grid=new_grid)
                    ds['sithick'] = ds['sithick'].where(lambda x:np.abs(x)<10)
                    sia_grid = self._weighted_sum(ds=ds['siconc'],weight=area,vid='siconc',calc=False)
                    ds_sum = self._weighted_sum(ds=ds[vids['aswgt'][0]],weight=sia_grid,vid=vids['aswgt'][0])
                    ds_sum = self._remove_vars(ds_sum)
                    datasets.append(ds_sum)

                if sid in source_ids_cloud:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat['aswgt'].search(source_id=sid)
                        if self.members=='first':
                            if subset.search(member_id='r1i1p1f1').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f1')
                            elif subset.search(member_id='r1i1p1f2').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f2')
                            elif subset.search(member_id='r1i1p1f3').df.member_id.size>=1:
                                subset = subset.search(member_id='r1i1p1f3')
                            else:
                                subset = subset.search(member_id='r1i1p2f1')
                            agg=True
                        ds_dict = subset.to_dataset_dict(aggregate=agg
                                                         ,xarray_open_kwargs={"consolidated": True,'decode_times': True,'use_cftime':True,'chunks': {'time': self.chunks}}
                                                                    ,storage_options={"anon": True},preprocess=complete_preprocessing,progressbar=False
                                                        ,xarray_combine_by_coords_kwargs=xr_dict)
                        dsets = list(map(self._postprocessing,ds_dict.values()))
                        dsets = list(map(self._remove_vars,dsets))
                        #dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                        #    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #    for exp in self.experiment_id]
                        #dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                        #    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                        #    for exp in self.experiment_id]
                        dsets = [xr.concat([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                            ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override',dim='member_id',join='outer') 
                            for exp in self.experiment_id]
                        [add_aswgt(ds) for ds in dsets if len(ds.data_vars)!=0 and 'x' and 'y' in ds.coords]
                        
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)                 
                            print('Warning: {} from {} did not load (grid={})'.format(vids['aswgt'],sid,grid))
                        continue

                if sid in source_ids_esgf:
                    if sid in self.skip_sids:
                        continue
                    try:
                        subset = cat_esgf['aswgt'].clone().search(**ChainMap({'source_id':sid}, cat_esgf['awgt'].last_search),quiet=True
                                                               ).remove_incomplete(complete=lambda sub_df: filter_missing(sub_df, missing=cat_esgf['aswgt'].df))
                        if self.members=='first':
                            subset.remove_ensembles()
                    
                        #ds_dict = subset.to_dataset_dict(prefer_streaming=True, add_measures=False, quiet=True)
                        #dsets = list(map(complete_preprocessing,ds_dict.values()))
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = load_from_catalog(
                            catalog=subset,
                            chunks = {'time':self.chunks},
                            preprocess=complete_preprocessing,
                            postprocess=self._postprocessing,
                            prefer_opendap=False,
                            combine_method='manual',
                            esgf_url = self.esgf_url
                        )
                        #dsets = list(map(self._postprocessing,dsets))
                        dsets = list(map(self._remove_vars,dsets))
                        try:
                            dsets = [xr.combine_by_coords([ds for ds in dsets if ds.experiment_id.values[0]==exp]
                                    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                                    for exp in self.experiment_id]
                        except Exception as error:
                            dsets = [xr.combine_nested([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                                    ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                                    for exp in self.experiment_id]
                        [add_aswgt(ds) for ds in dsets if len(ds.data_vars)!=0 and 'x' and 'y' in ds.coords]
                        
                    except Exception as error:
                        if self.verbose==True:
                            print("An error occurred:", type(error).__name__, "–", error)
                            print('Warning:{} from {} did not load (grid={})'.format(vids['aswgt'],sid,grid))
                        continue   
                try:
                    if sid in unstr_sids:
                        new_grid = grid_CESM2
                        area = grid_CESM2.areacello 
                    else:
                        area = self._get_area(sid,grid)
                        try:
                            area.compute()
                        except Exception:
                            area = calc_areacello(dsets[0])
                        if area is None or sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1','NESM3']:
                            #if sid in ['KIOST-ESM','CAS-ESM2-0','BCC-CSM2-MR','BCC-ESM1']:
                            #    dset,area = self._fix_grids(dsets[0],area.areacello)
                            area = calc_areacello(dsets[0])
                        area = area.areacello
                except (IndexError,AttributeError,ValueError):
                    if self.verbose==True:
                        print('Warning: area from {} not found (grid={})'.format(sid,grid))
                    continue
                    
        try:
            dsets = list(map(self._convert_to_dataset,datasets))
            ds = xr.concat([xr.concat([ds for ds in dsets if ds.experiment_id.values[0]==exp],'member_id'
                            ,coords='minimal',data_vars='minimal',compat='override',combine_attrs='override') 
                       for exp in self.experiment_id],'experiment_id',combine_attrs='override')
            return ds
            #return dsets
        except ValueError:
            pass

    def load_data(self):
        if self.source_id is None:
            if self.verbose==True:
                print('Warning: No source_id given')
            return 
        vids,tid=self._get_vids(),self.tid
        datasets = []
        for grid in (pbar := tqdm(self.grid_label)):
            pbar.set_description(f"Processing data for each grid")
            cat,cat_esgf = self._get_cat(vids=vids,tid=tid,grid=grid)
            if self.client==True:
                data = delayed(self._load_data)(vids=vids,cat=cat,cat_esgf=cat_esgf,tid=tid,grid=grid)
                data = data.compute()
            else:
                data = self._load_data(vids=vids,cat=cat,cat_esgf=cat_esgf,tid=tid,grid=grid)
            if data is None:
                continue
            #if len(data.data_vars) != 0:
            #    datasets.append(data)
            if len(data) != 0:
                datasets.append(data)
        if len(datasets) == 0:
            if self.verbose==True:
                print('Warning: No data found for grid label: {}'.format(grid))
            return
        else:
            return xr.concat(datasets,'member_id').drop_duplicates('member_id')
            #return datasets
    
    def data_summary(self):
        vids,tid=self._get_vids(),self.tid
        cat,cat_esgf = self._get_cat(vids=vids,tid=tid,grid=self.grid_label[0])
        return cat,cat_esgf


def shift_and_convert_time(ds, delta_years=1700):
    # Calculate year difference from shift point
    #original_years = xr.DataArray([t.year for t in ds.time.values], dims="time")
    #min_year = int(original_years.min())
        
    # If already >= shift_years, no need to shift
    #if min_year >= shift_years:
    #    return ds.convert_calendar("proleptic_gregorian", use_cftime=False)
        
    # Shift the year values
    #delta_years = shift_years - min_year
    shifted_times = [pd.Timestamp(f"{t.year + delta_years}-{t.month:02d}-{t.day:02d}") for t in ds.time.values]
        
    # Replace time coordinate
    ds = ds.assign_coords(time=("time", shifted_times))
        
    # Optionally save metadata about the shift
    ds.attrs["time_shift_years"] = delta_years
    return ds

class CESM():
    def __init__(self,variable, source_id='CESM2-LENS', sector_mean=None,sector_sum=None,time_chunks=200, shift_years=None,shift_months=-1,
                 new_grid=None, method='conservative_normed',sic_mask=None,verbose=False,data_path='/',
                 opendap_forcing='cmip6', opendap_experiment=None, opendap_time_range=None):
        #at the moment, there is no ability to take spatial averages and sums here, but maybe there will be a need to add in the future
        self.variable = variable
        self.source_id = source_id
        self.time_chunks = time_chunks
        self.data_path = data_path
        self.new_grid = new_grid
        self.method = method
        self.data_path=data_path
        self.sic_mask = sic_mask
        self.verbose = verbose
        self.shift_years = shift_years
        self.shift_months = shift_months
        self.opendap_forcing = opendap_forcing      # 'cmip6' or 'smbb'
        self.opendap_experiment = opendap_experiment  # 'BHIST', 'BSSP370', or None (both)
        self.opendap_time_range = opendap_time_range  # (start_year, end_year) or None

    def preprocess_static_vars(self,ds):
        for var in list(ds.data_vars):
            if var not in ['hs','hi','sifb','aice','sisnthick','sithick','apeff_ai','apond_ai','siflswutop','siflswdtop','siflswdbot','uvel_d','vvel_d',
                          'TREFHT','sitemptop', "sidmassgrowthbot", "sidmassgrowthwat", "sidmasssi", "sidmassmelttop", "sidmassmeltbot", "sidmasslat",
                           "sidmassevapsubl","sidmassdyn","sndmassmelt","sndmasssnf","sndmassubl",'sifllwdtop','sidmassth']:
                ds = ds.set_coords(var)
        return ds

    def rename(self,ds):
        if 'TLON' in ds.coords:
            ds = ds.rename({'TLON':'lon','TLAT':'lat'})
        if 'nj' in ds.coords:
            ds = ds.rename({'ni':'x','nj':'y'})
        if 'nvertices' in ds.coords:
            ds=ds.rename({'nvertices':'vertex'})
        if 'vertex' in ds.variables and 'y' in ds.variables:
            if 'lont_bounds' in ds.variables:
                lon_corners = cfxr.bounds_to_vertices(ds.lont_bounds.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
                lat_corners = cfxr.bounds_to_vertices(ds.latt_bounds.chunk(dict(vertex=-1,y=-1,x=-1)), "vertex", order=None)
            ds=ds.assign_coords(lon_b=lon_corners, lat_b=lat_corners)
            ds=ds.rename({'y_vertices':'y_b','x_vertices':'x_b'})
        if 'sisnthick' in ds.variables:
            ds=ds.rename({'sisnthick':'snt'})
        if 'aice' in ds.variables:
            ds=ds.rename({'aice':'sic'})
        if 'sithick' in ds.variables:
            ds=ds.rename({'sithick':'sit'})
        if 'TREFHT' in ds.variables:
            ds=ds.rename({'TREFHT':'tas'})
        #if 'apond_ai' in ds.variables:
        #    ds=ds.rename({'apond_ai':'simpconc'})
        if 'tmask' in ds.coords:
            ds=ds.rename({'tmask':'mask'})
        if 'tarea' in ds.coords:
            ds=ds.rename({'tarea':'areacello'})
        return ds
        
    def _extract_member_number(self,casename):
        # Try LE2-style:
        match_le2 = re.search(r"LE2-(\d+)\.(\d+)", casename)
        if match_le2:
            return match_le2.group(1) + match_le2.group(2)
    
        # Try CMIP6-style:
        match_cmip6 = re.search(r"CMIP6-[^.]*\.\d+\.(\d+)\.cice", casename)
        if match_cmip6:
            return match_cmip6.group(1)

        return None
    
    def _extract_experiment_id(self,casename):
        parts = casename.split(".")
        if len(parts) >= 3:
            return parts[2]  # 0=b, 1=e21, 2=BSSP370cmip6
        return None
    
    def _standardize_experiment_id(self,casename):
        exp = self._extract_experiment_id(casename)
        if exp:
            match = re.search(r"(SSP\d{3})", exp.upper())
            return match.group(1).lower() if match else exp.lower()
        return None
    
    def update_member_id(self,ds):
        ds = ds.expand_dims({'member_id':[self.source_id+'_'+self._extract_member_number(ds.attrs['title'])]})    
        ds['member_id'] = ds.member_id.astype('<U25')
        ds = ds.expand_dims({'experiment_id':[self._standardize_experiment_id(ds.attrs['title'])]})
        ds['experiment_id'] = ds.experiment_id.astype('<U10')
        return ds

    def update_member_id(self, ds):
        # Use the real filename from encoding
        filename = os.path.basename(ds.encoding.get('source', ''))
        #print(f"DEBUG: using filename = {filename}")
    
        member = self._extract_member_number(filename)
        if not member:
            ds = ds.expand_dims({'member_id': [f"{self.source_id}"]})
            ds['member_id'] = ds.member_id.astype('<U25')

        else:
            ds = ds.expand_dims({'member_id': [f"{self.source_id}_{member}"]})
            ds['member_id'] = ds.member_id.astype('<U25')
            
        experiment = self._standardize_experiment_id(filename)
        if not experiment:
            raise ValueError(f"Could not extract experiment ID from: {filename}")
    
        ds = ds.expand_dims({'experiment_id': [experiment]})
        ds['experiment_id'] = ds.experiment_id.astype('<U10')
    
        return ds
    
    def convert_time(self,ds):
        try:
            exp_id = ds.experiment_id
            # Apply time slicing only if it makes sense
            if exp_id == "historical":
                ds = ds.sel(time=slice("1850-01", "2014-12"))
            elif exp_id in ['ssp126','ssp245','ssp370','ssp585']:
                ds = ds.sel(time=slice("2015-01", "2100-12"))

            if self.shift_years is not None:
                ds = shift_and_convert_time(ds,self.shift_years)
    
            # Fix time formatting
            if "time" in ds.coords:
                ds["time"] = pd.to_datetime(ds.time.dt.strftime("%Y-%m-%d").values)
    
            # Drop duplicate timestamps
            ds = ds.drop_duplicates("time", keep="first")
            return ds
    
        except Exception as e:
            print(f"[convert_time] ⚠️ Error: {e}")
            return ds
    
    def shift_time_by_months(self,ds, n=1):
        from dateutil.relativedelta import relativedelta
        """
        Shift the time coordinate of an xarray Dataset by `n` months.
        Handles both cftime and datetime64 time types.
        """
        try:
            time = ds["time"].values
            first_time = time[0]
    
            # Use cftime logic
            if isinstance(first_time, cftime.DatetimeNoLeap) or isinstance(first_time, cftime.Datetime360Day) or isinstance(first_time, cftime.DatetimeJulian):
                shifted_time = [t + relativedelta(months=n) for t in time]
    
            # Use pandas datetime logic
            elif isinstance(first_time, (np.datetime64, pd.Timestamp)):
                shifted_time = pd.to_datetime(time) + pd.DateOffset(months=n)
    
            else:
                raise TypeError(f"Unsupported time type: {type(first_time)}")
    
            # Assign back to dataset
            ds = ds.copy()
            ds["time"] = ("time", shifted_time)
            return ds
    
        except Exception as e:
            print(f"[shift_time_by_months] ⚠️ Error: {e}")
            return ds

    def drop_aux_lonlat(self, ds):
        return ds.drop_vars(["zlon", "zlat", "zlon_bnds", "zlon_bnds"], errors="ignore")

    # ── OPeNDAP / THREDDS helpers ──────────────────────────────────────────────

    def _data_path_valid(self) -> bool:
        """Return True if data_path (used as a glob prefix) matches files for self.variable."""
        if not self.data_path:
            return False
        return bool(glob.glob(f"{self.data_path}*{self.variable}*.nc"))

    def _in_time_range(self, filename: str) -> bool:
        """Return True if the file's time range overlaps with opendap_time_range."""
        if self.opendap_time_range is None:
            return True
        m = re.search(r'\.(\d{6})-(\d{6})\.nc$', filename)
        if not m:
            return True  # can't determine, include it
        file_start = int(m.group(1))   # YYYYMM
        file_end   = int(m.group(2))   # YYYYMM
        req_start  = int(self.opendap_time_range[0]) * 100 + 1
        req_end    = int(self.opendap_time_range[1]) * 100 + 12
        return file_start <= req_end and file_end >= req_start

    def _list_opendap_urls(self, variable: str = None) -> list:
        """
        List OPeNDAP URLs for `variable` by scraping the THREDDS HTML catalog.

        Scraping HTML is more reliable than parsing the catalog XML, which uses
        a DatasetScan stub that doesn't embed child datasets in the XML response.

        Files are stored under month_1/{variable}/ — e.g.
        month_1/sidmassth/b.e21.BHISTcmip6.f09_g17.LE2-1011.001.cice.h.sidmassth.185001-185912.nc
        """
        var = variable if variable is not None else self.variable
        forcing = self.opendap_forcing.lower()
        exclude = 'smbb' if forcing == 'cmip6' else 'cmip6'
        experiment = self.opendap_experiment  # 'BHIST', 'BSSP370', or None

        var_path = f"{_CESM2LE_THREDDS_PATH}/{var}"
        catalog_url = f"{_CESM2LE_THREDDS_BASE}/catalog{var_path}/catalog.html"

        if catalog_url not in _THREDDS_CATALOG_CACHE:
            last_err = None
            for attempt in range(1, 5):
                try:
                    resp = requests.get(catalog_url, timeout=120)
                    resp.raise_for_status()
                    _THREDDS_CATALOG_CACHE[catalog_url] = resp.text
                    break
                except Exception as e:
                    last_err = e
                    if attempt < 4:
                        time.sleep(2 ** attempt)
            else:
                raise last_err
        html = _THREDDS_CATALOG_CACHE[catalog_url]

        # Extract all .nc filenames from the HTML page
        all_nc = sorted(set(re.findall(r'[^\s"<>/]+\.nc', html)))

        urls = []
        for name in all_nc:
            if (var in name
                    and forcing in name
                    and exclude not in name
                    and (experiment is None or experiment in name)
                    and self._in_time_range(name)):
                urls.append(f"{_CESM2LE_THREDDS_BASE}/dodsC{var_path}/{name}")

        if not urls and self.verbose:
            has_var     = [n for n in all_nc if var in n]
            has_forcing = [n for n in has_var if forcing in n and exclude not in n]
            print(f"  [OPeNDAP] catalog (HTML): {catalog_url}")
            print(f"  [OPeNDAP] .nc files found in page   : {len(all_nc)}")
            print(f"  [OPeNDAP] contain '{var}'            : {len(has_var)}")
            print(f"  [OPeNDAP] match forcing '{forcing}'  : {len(has_forcing)}")
            print(f"  [OPeNDAP] match experiment '{experiment}': "
                  f"{len([n for n in has_forcing if experiment is None or experiment in n])}")
            if has_var and not has_forcing:
                labels = sorted({n.split('.')[2] for n in has_var if '.' in n})
                print(f"  [OPeNDAP] experiment labels found: {labels}")
                print(f"  [OPeNDAP] hint: try opendap_forcing='smbb'")
            if has_forcing and experiment:
                labels = sorted({n.split('.')[2] for n in has_forcing if '.' in n})
                print(f"  [OPeNDAP] experiment labels after forcing filter: {labels}")
                print(f"  [OPeNDAP] hint: try opendap_experiment='BHIST' or 'BSSP370'")

        return sorted(urls)

    def inspect_catalog(self, variable: str = None) -> None:
        """Print filenames found in the THREDDS HTML catalog for debugging."""
        var = variable if variable is not None else self.variable
        url = (f"{_CESM2LE_THREDDS_BASE}/catalog"
               f"{_CESM2LE_THREDDS_PATH}/{var}/catalog.html")
        resp = requests.get(url, timeout=60)
        print(f"Status: {resp.status_code}  URL: {url}\n")
        names = sorted(set(re.findall(r'[^\s"<>/]+\.nc', resp.text)))
        print(f"Found {len(names)} .nc filenames:")
        for n in names[:20]:
            print(f"  {n}")
        if len(names) > 20:
            print(f"  ... ({len(names) - 20} more)")

    def _open_via_opendap(self, variable: str = None) -> xr.Dataset:
        """Open CESM2-LE files for `variable` via OPeNDAP using the THREDDS catalog."""
        var = variable if variable is not None else self.variable
        # Force verbose diagnostics on this call so the filter counts print on failure
        _prev_verbose, self.verbose = self.verbose, True
        urls = self._list_opendap_urls(variable=var)
        self.verbose = _prev_verbose
        if not urls:
            raise FileNotFoundError(
                f"No CESM2-LE OPeNDAP files found for variable '{var}' "
                f"with forcing='{self.opendap_forcing}'. "
                f"See diagnostic output above."
            )
        if self.verbose:
            print(f"  OPeNDAP: opening {len(urls)} files sequentially for '{var}'")
        xr_kw = {
            'coords': 'minimal', 'data_vars': 'minimal',
            'compat': 'override', 'combine_attrs': 'override',
        }

        datasets = []
        for url in (pbar (pbar := tqdm(urls, desc=f"OPeNDAP {var}", unit="file", leave=False)):
            pbar.set_postfix_str(os.path.basename(url)[:40])= tqdm(urls, desc=f"OPeNDAP {var}", unit="file", leave=False)):
            pbar.set_postfix_str(os.path.basename(url)[:40])
            for attempt in range(1, 5):
                try:
                    ds = xr.open_dataset(
                        url,
                        chunks={'time': self.time_chunks},
                        decode_cf=True,
                        use_cftime=True,
                        engine='netcdf4',
                    )
                    ds = self.preprocessing(ds)
                    datasets.append(ds)
                    break
                except Exception as e:
                    if attempt < 4:
                        if self.verbose:
                            print(f"  Retry {attempt}/4 for {os.path.basename(url)}: "
                                  f"{type(e).__name__}")
                        time.sleep(2 ** attempt)
                    elif self.verbose:
                        print(f"  Warning: skipping {os.path.basename(url)} "
                              f"after 4 attempts: {type(e).__name__}: {e}")

        if not datasets:
            raise OSError(f"No files could be opened for variable '{var}'")

        return xr.combine_by_coords(datasets, **xr_kw)

    # ── end OPeNDAP helpers ────────────────────────────────────────────────────

    def preprocessing(self,ds):
        ds = self.preprocess_static_vars(ds)
        ds = self.drop_aux_lonlat(ds)
        ds = promote_empty_dims(ds)
        ds = correct_coordinates(ds)
        ds = self.rename(ds)
        ds = self.update_member_id(ds)
        ds = self.convert_time(ds)
        ds = self.shift_time_by_months(ds, n=self.shift_months)
        return ds
    
    def calc_freeboard(self,ds):
        rho_w = 1026 #kg/m3
        rho_i = 916 #kg/m3
        rho_sn = 330 #kg/m3
        H_i = ds.sit.where(lambda x:np.abs(x)<10)
        H_sn = ds.snt.where(lambda x:np.abs(x)<10)
        ds['sifb_d'] = H_i * ((rho_w - rho_i)/rho_w) - H_sn * (rho_sn/rho_w)
        return ds

    def calc_sit(self,ds):
        ds['sit'] = ds.hi/ds.sic
        return ds

    def calc_snt(self,ds):
        ds['snt'] = ds.hs/ds.sic
        return ds
    
    def calc_simpconc(self,ds):
        ds['simpconc'] = ds[self.variable]*100
        return ds

    def regrid(self,ds,new_grid=None):
        if self.new_grid == None:
            ds_out=new_grid
        else:
            ds_out=self.new_grid
        locstream_in=False
        skipna=True
        extrap_method=None
        regridder = xe.Regridder(ds, ds_out, self.method
                                ,ignore_degenerate=True,extrap_method=extrap_method,periodic=True,locstream_in=locstream_in)
        ds_new = regridder(ds,skipna=skipna)
        #return ds_new.where(ds_new.mask).chunk({'y':ds_new.y.size,'x':ds_new.x.size})
        if 'mask' in ds_new.coords:
            return ds_new.where(ds_new.mask)
        else:
            return ds_new
    
    def _open_vars(self, vars_to_load):
        paths = []
        for var in vars_to_load:
            paths.extend(glob.glob(f"{self.data_path}*{var}*.nc"))
    
        if not paths:
            raise FileNotFoundError(f"No files found for variables: {vars_to_load}")
    
        xr_dict = {
            'coords': 'minimal',
            'data_vars': 'minimal',
            'compat': 'override',
            'combine_attrs': 'override',
        }
    
        ds = xr.open_mfdataset(
            sorted(paths),
            chunks={'time': self.time_chunks},
            decode_cf=True,
            use_cftime=True,
            preprocess=self.preprocessing,
            **xr_dict
        )
        return ds    

    def load_data(self):
        xr_dict = {'coords':'minimal','data_vars':'minimal','compat':'override','combine_attrs':'override'}

        # Fall back to THREDDS OPeNDAP when the local data_path is unavailable
        use_opendap = (self.source_id == 'CESM2-LE' and not self._data_path_valid())
        if use_opendap and self.verbose:
            print(f"data_path '{self.data_path}' not valid; falling back to THREDDS OPeNDAP")

        if self.variable in ['sifb','aice','sisnthick','sithick','siflswutop','siflswdtop','siflswdbot','uvel_d','vvel_d',
                          'TREFHT','sitemptop', "sidmassgrowthbot", "sidmassgrowthwat", "sidmasssi", "sidmassmelttop", "sidmassmeltbot", "sidmasslat",
                           "sidmassevapsubl","sidmassdyn","sndmassmelt","sndmasssnf","sndmassubl",'sifllwdtop','sidmassth']:
            if use_opendap:
                ds = self._open_via_opendap()
            else:
                ds = xr.open_mfdataset(self.data_path+'*{}*.nc'.format(self.variable)
                                   ,chunks={'time':self.time_chunks}
                                , decode_cf=True, use_cftime=True,preprocess=self.preprocessing,**xr_dict)
            #ds = ds.chunk({'time': self.time_chunks})
        elif self.variable in ['hi','hs']:
            if use_opendap:
                ds = self._open_via_opendap()
            else:
                ds = xr.open_mfdataset(self.data_path+'*.nc'
                                   ,chunks={'time':self.time_chunks}
                                , decode_cf=True, use_cftime=True,preprocess=self.preprocessing,**xr_dict)
            #ds = ds.chunk({'time': self.time_chunks})
            ds = self.calc_sit(ds)
            ds = self.calc_snt(ds)
        elif self.variable == 'sifb_d':
            if use_opendap:
                ds = xr.merge([
                    self._open_via_opendap('sithick'),
                    self._open_via_opendap('sisnthick'),
                    self._open_via_opendap('aice'),
                ])
            else:
                ds = self._open_vars(['sithick', 'sisnthick', 'aice'])
            if 'hi' in ds.data_vars: #not needed unless sithick and sisnthick are missing, then use hi/hs and convert
                ds = self.calc_sit(ds)
                ds = self.calc_snt(ds)
            ds = self.calc_freeboard(ds)
        elif self.variable in ['apond_ai','apeff_ai']:
            if use_opendap:
                ds = self._open_via_opendap()
            else:
                ds = xr.open_mfdataset(self.data_path+'*{}*.nc'.format(self.variable)
                                   ,chunks={'time':self.time_chunks}
                                , decode_cf=True, use_cftime=True,preprocess=self.preprocessing,**xr_dict)
            #ds = ds.chunk({'time': self.time_chunks})
            ds = self.calc_simpconc(ds)
        else:
            print('variable not available')
            return
        if self.sic_mask!=None:
            if 'sic' not in ds.data_vars:
                if use_opendap:
                    sic = self._open_via_opendap('aice').sic
                    sic = sic.chunk({'time': self.time_chunks})
                else:
                    sic = xr.open_mfdataset(self.data_path+'*{}*.nc'.format('aice')
                                   ,chunks={'time':1}
                                , decode_cf=True, use_cftime=True,preprocess=self.preprocessing,**xr_dict).sic
                    sic = sic.chunk({'time': self.time_chunks})
            else:
                sic = ds.sic
            ds = ds.where(sic>self.sic_mask)

        if self.variable=='sithick':
            self.variable = 'sit'
        if self.variable=='sisnthick':
            self.variable = 'snt'
        if self.variable=='aice':
            self.variable = 'sic'
        if self.variable in ['apond_ai','apeff_ai']:
            self.variable = 'simpconc'
        if self.variable=='hi':
            self.variable = 'sit'
        if self.variable=='hs':
            self.variable = 'snt'
        if self.variable=='TREFHT':
            self.variable = 'tas'
        
            
        keep_coords = ['experiment_id','member_id','lon','lat','x','y','mask','areacello','time']

        if self.new_grid==None:
            ds = ds.drop_vars([var for var in ds.data_vars if var != self.variable])
            return ds.drop_vars([coord for coord in ds.coords if coord not in keep_coords])
        else:
            ds_regridded = self.regrid(ds.drop_vars([var for var in ds.data_vars if var != self.variable]))
            return ds_regridded.drop_vars([coord for coord in ds_regridded.coords if coord not in keep_coords])